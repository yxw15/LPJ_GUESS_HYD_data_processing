# ==========================================================================
# CROSS-VALIDATION: Gc & Leaf Water Potential (ΨL)
# Validates that parameter selection is robust, not just overfitting.
#
# Strategies:
#   For Gc (daily, many obs):
#     1. Leave-One-Month-Out (LOMO): 4-fold, train on 3 months, test on 1
#     2. Temporal 80/20 block split: train on early season, test on late
#     3. Leave-One-Treatment-Out (LOTO): train on control ↔ test on drought
#
#   For ΨL (~6 obs per species×treatment):
#     4. Leave-One-Out (LOO): train on 5 points, test on 1
#     5. Pooled-species LOO: train on all species except one date
#
# Key question: Do the "best" parameters identified on training data
#               still perform well on held-out test data?
#   YES → model is robust, parameters are identifiable
#   NO  → overfitting, parameter estimates are not reliable
#
# Usage:
#   Rscript cross_validation_Gc_PSI_hydraulic.R
#
# Output:
#   Figures/lpj_guess_stem_storage/cross_validation/
#     cv_summary_Gc.csv              — per-fold train/test metrics
#     cv_summary_PSI.csv             — per-fold train/test metrics
#     cv_stability.csv               — how often same params are selected
#     cv_Gc_scatter.png              — train vs test R² / KGE
#     cv_PSI_scatter.png             — train vs test R² / KGE
#     cv_parameter_stability.png     — which params are consistently selected
# ==========================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(purrr)

# ==========================================================================
# 0. PATHS & SETUP
# ==========================================================================
BASE_DIR       <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
RESULTS_BASE   <- file.path(BASE_DIR, "results_lpj/results_sensitivity")
OUTPUT_DIR     <- file.path(BASE_DIR, "Figures/lpj_guess_stem_storage/cross_validation")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(BASE_DIR)

# ==========================================================================
# 1. METRIC FUNCTIONS (same as sensitivity analysis)
# ==========================================================================
safe_cor <- function(x, y) {
  if (length(x) < 3 || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0)
    return(NA_real_)
  cor(x, y, use = "complete.obs")
}

compute_kge <- function(sim, obs) {
  valid <- complete.cases(sim, obs)
  if (sum(valid) < 3) return(NA_real_)
  sim_v <- sim[valid]
  obs_v <- obs[valid]
  if (sd(obs_v) == 0 || mean(obs_v) == 0) return(NA_real_)

  r     <- cor(sim_v, obs_v)
  alpha <- sd(sim_v) / sd(obs_v)
  beta  <- mean(sim_v) / mean(obs_v)
  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

compute_metrics <- function(sim, obs) {
  valid <- complete.cases(sim, obs)
  if (sum(valid) < 3) {
    return(data.frame(
      n = sum(valid), pearson_r = NA_real_, r_squared = NA_real_,
      rmse = NA_real_, nrmse_pct = NA_real_, bias = NA_real_,
      slope = NA_real_, kge = NA_real_
    ))
  }
  s <- sim[valid]
  o <- obs[valid]

  r_val   <- cor(s, o)
  fit     <- lm(s ~ o)
  sl      <- coef(fit)[2]
  rmse_v  <- sqrt(mean((s - o)^2))
  nrmse_v <- (rmse_v / mean(o)) * 100
  bias_v  <- mean(s - o)
  kge_v   <- compute_kge(s, o)

  data.frame(
    n          = sum(valid),
    pearson_r  = r_val,
    r_squared  = r_val^2,
    rmse       = rmse_v,
    nrmse_pct  = nrmse_v,
    bias       = bias_v,
    slope      = sl,
    kge        = kge_v
  )
}

# Composite score (same weighting as sensitivity script)
compute_composite <- function(r_squared, kge, slope, bias, nrmse_pct) {
  0.35 * r_squared +
    0.30 * kge +
    0.15 * (1 - pmin(abs(1 - slope), 1)) +
    0.10 * (1 - pmin(abs(bias), 1)) +
    0.10 * (1 - pmin(nrmse_pct / 200, 1))
}

# ==========================================================================
# 2. READ LPJ OUTPUT (same as sensitivity script)
# ==========================================================================
read_lpj_run <- function(run_dir, species_pft) {
  dgc_file      <- file.path(run_dir, "dgc.out")
  dpsileaf_file <- file.path(run_dir, "dpsileaf.out")

  if (!file.exists(dgc_file) || !file.exists(dpsileaf_file)) return(NULL)

  gc_data <- tryCatch(read.table(dgc_file, header = TRUE, check.names = FALSE),
                      error = function(e) NULL)
  pl_data <- tryCatch(read.table(dpsileaf_file, header = TRUE, check.names = FALSE),
                      error = function(e) NULL)

  if (is.null(gc_data) || is.null(pl_data)) return(NULL)

  colnames(gc_data) <- trimws(colnames(gc_data))
  colnames(pl_data) <- trimws(colnames(pl_data))

  if (!species_pft %in% colnames(gc_data) || !species_pft %in% colnames(pl_data))
    return(NULL)

  pft_name <- unname(species_pft)
  colnames(gc_data)[colnames(gc_data) == pft_name] <- "Gc"
  gc_data <- gc_data %>%
    mutate(date = as.Date(Day, origin = paste0(Year, "-01-01"))) %>%
    select(date, Gc)

  colnames(pl_data)[colnames(pl_data) == pft_name] <- "psi_leaf"
  pl_data <- pl_data %>%
    mutate(date = as.Date(Day, origin = paste0(Year, "-01-01"))) %>%
    select(date, psi_leaf)

  inner_join(gc_data, pl_data, by = "date") %>%
    filter(date >= as.Date("2023-01-01"))
}

# ==========================================================================
# 3. READ OBSERVATIONS
# ==========================================================================
obs_gc <- tryCatch({
  read.csv(file.path(BASE_DIR, "SCCII/sap_flux_gc_daytime_Climate_filter.csv")) %>%
    mutate(
      date      = as.Date(date),
      species   = tolower(species),
      treatment = case_when(treatment == "control" ~ "control",
                            treatment == "treatment" ~ "drought")
    ) %>%
    filter(date >= as.Date("2023-01-01"), !is.na(G_ms), G_ms <= 12)
}, error = function(e) NULL)

obs_psi <- tryCatch({
  bind_rows(
    read.csv(file.path(BASE_DIR, "SCCII/psiL_hoelstein_drought.csv")) %>%
      mutate(treatment = "drought"),
    read.csv(file.path(BASE_DIR, "SCCII/psiL_hoelstein_control.csv")) %>%
      mutate(treatment = "control")
  ) %>%
    mutate(date = as.Date(date), month = month(date)) %>%
    filter(month >= 6, month <= 9) %>%
    rename(species = species_name) %>%
    group_by(date, species, treatment) %>%
    summarise(psi_leaf_md_obs = mean(md_wp_av, na.rm = TRUE), .groups = "drop")
}, error = function(e) NULL)

# ==========================================================================
# 4. MAPPINGS & PARAMETER DEFINITIONS
# ==========================================================================
species_pft_map <- c(
  Oak    = "Que_rob",
  Beech  = "Fag_syl",
  Spruce = "Pic_abi",
  Pine   = "Pin_syl"
)
species_lower_map <- c(
  Oak    = "oak",
  Beech  = "beech",
  Spruce = "spruce",
  Pine   = "pine"
)

parameter_bases <- list(
  Oak = list(
    control = list(kr = 1.8, ks = 2.1, kl = 10,   iso = 0.4, vc = 200, theta = 0.55),
    drought = list(kr = 1.8, ks = 2.1, kl = 10,   iso = 0.4, vc = 150, theta = 0.65)
  ),
  Beech = list(
    control = list(kr = 0.8, ks = 1.0, kl = 25,   iso = 0.6, vc = 200, theta = 0.55),
    drought = list(kr = 0.8, ks = 1.0, kl = 25,   iso = 0.6, vc = 150, theta = 0.70)
  ),
  Spruce = list(
    control = list(kr = 0.5, ks = 0.7, kl = 33.1, iso = 0.7, vc = 300, theta = 0.55),
    drought = list(kr = 0.5, ks = 0.7, kl = 33.1, iso = 0.7, vc = 300, theta = 0.55)
  ),
  Pine = list(
    control = list(kr = 0.6, ks = 0.8, kl = 12.5, iso = 0.7, vc = 300, theta = 0.60),
    drought = list(kr = 0.6, ks = 0.8, kl = 12.5, iso = 0.7, vc = 300, theta = 0.60)
  )
)

oat_kr    <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5)
oat_ks    <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9)
oat_kl    <- c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35)
oat_iso   <- c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
oat_vc    <- c(150, 200, 250, 300, 350, 400)
oat_theta <- c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70)

# Build parameter combo list (same as sensitivity script)
build_combos <- function(bases) {
  combos <- list()
  fmt_tag_1dp <- function(x) {
    x <- round(x, digits = 8)
    ifelse(x == round(x), sprintf("%.1f", x), as.character(x))
  }
  fmt_tag_2dp <- function(x) sprintf("%.2f", round(x, digits = 8))

  # Baseline
  combos[[1]] <- list(
    tag = "baseline",
    kr_max = bases$kr, ks_max = bases$ks,
    kl_max = bases$kl, isohydricity = bases$iso,
    volumetric_capacitance = bases$vc, sapwood_theta_sat = bases$theta
  )
  idx <- 2

  # kr_max OAT
  for (v in oat_kr) {
    combos[[idx]] <- list(tag = paste0("kr", fmt_tag_1dp(v)),
      kr_max = v, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta)
    idx <- idx + 1
  }
  # ks_max OAT
  for (v in oat_ks) {
    combos[[idx]] <- list(tag = paste0("ks", fmt_tag_1dp(v)),
      kr_max = bases$kr, ks_max = v, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta)
    idx <- idx + 1
  }
  # kl_max OAT
  for (v in oat_kl) {
    combos[[idx]] <- list(tag = paste0("kl", as.character(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = v,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta)
    idx <- idx + 1
  }
  # isohydricity OAT
  for (v in oat_iso) {
    combos[[idx]] <- list(tag = paste0("iso", fmt_tag_1dp(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = v, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta)
    idx <- idx + 1
  }
  # volumetric_capacitance OAT
  for (v in oat_vc) {
    combos[[idx]] <- list(tag = paste0("vc", as.character(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = v,
      sapwood_theta_sat = bases$theta)
    idx <- idx + 1
  }
  # sapwood_theta_sat OAT
  for (v in oat_theta) {
    combos[[idx]] <- list(tag = paste0("theta", fmt_tag_2dp(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = v)
    idx <- idx + 1
  }
  combos
}

# ==========================================================================
# 5. PRELOAD ALL LPJ MODEL OUTPUT (cache to avoid re-reading)
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  PRELOADING LPJ-GUESS OUTPUT\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Structure: lpj_cache[[sp]][[cond]][[tag]] = data.frame(date, Gc, psi_leaf)
lpj_cache <- list()

for (sp in names(parameter_bases)) {
  lpj_cache[[sp]] <- list()
  pft <- species_pft_map[sp]

  for (cond in c("control", "drought")) {
    lpj_cache[[sp]][[cond]] <- list()
    bases <- parameter_bases[[sp]][[cond]]
    combos <- build_combos(bases)

    for (cmb in combos) {
      dir_name <- sprintf("%s_%s_%s", sp, cond, cmb$tag)
      run_dir  <- file.path(RESULTS_BASE, cond, sp, dir_name)

      if (!dir.exists(run_dir)) {
        cat(sprintf("  [MISSING] %s\n", dir_name))
        next
      }

      lpj_data <- read_lpj_run(run_dir, pft)
      if (!is.null(lpj_data) && nrow(lpj_data) > 0) {
        lpj_cache[[sp]][[cond]][[cmb$tag]] <- lpj_data
      }
    }
    cat(sprintf("  %-8s %-10s: %d / %d runs loaded\n",
                sp, cond,
                length(lpj_cache[[sp]][[cond]]),
                length(combos)))
  }
}

# ==========================================================================
# 6. CROSS-VALIDATION: LEAVE-ONE-MONTH-OUT (LOMO) for Gc
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  CROSS-VALIDATION 1: LEAVE-ONE-MONTH-OUT (Gc)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Define growing season months to use as folds
growing_months <- 5:9  # May through September = 5 folds

cv_results_gc <- list()

for (sp in names(parameter_bases)) {
  sp_low <- species_lower_map[sp]

  for (cond in c("control", "drought")) {

    # Get observation dates for this species×condition
    obs_sub <- obs_gc %>%
      filter(species == sp_low, treatment == cond) %>%
      group_by(date) %>%
      summarise(G_obs = mean(G_ms, na.rm = TRUE), .groups = "drop") %>%
      mutate(month = month(date)) %>%
      filter(month %in% growing_months)

    if (nrow(obs_sub) < 20) {
      cat(sprintf("  [SKIP] %-8s %-10s: only %d obs days\n", sp, cond, nrow(obs_sub)))
      next
    }

    # For each test month
    for (test_month in growing_months) {

      train_dates <- obs_sub %>% filter(month != test_month) %>% pull(date)
      test_dates  <- obs_sub %>% filter(month == test_month) %>% pull(date)

      if (length(test_dates) < 5) next  # need minimum test points

      # Evaluate ALL parameter combos on TRAIN data
      best_train_score <- -Inf
      best_tag         <- NA_character_
      best_test_metrics <- NULL

      for (tag in names(lpj_cache[[sp]][[cond]])) {
        lpj_data <- lpj_cache[[sp]][[cond]][[tag]]
        if (is.null(lpj_data)) next

        # Training performance
        train_joined <- lpj_data %>%
          filter(date %in% train_dates) %>%
          inner_join(obs_sub %>% select(date, G_obs), by = "date") %>%
          filter(Gc <= 12)

        if (nrow(train_joined) < 5) next

        train_m <- compute_metrics(train_joined$Gc, train_joined$G_obs)

        train_score <- compute_composite(
          train_m$r_squared, train_m$kge,
          train_m$slope, train_m$bias, train_m$nrmse_pct
        )

        # Test performance
        test_joined <- lpj_data %>%
          filter(date %in% test_dates) %>%
          inner_join(obs_sub %>% select(date, G_obs), by = "date") %>%
          filter(Gc <= 12)

        test_m <- compute_metrics(test_joined$Gc, test_joined$G_obs)

        if (train_score > best_train_score) {
          best_train_score <- train_score
          best_tag         <- tag
          best_train_m     <- train_m
          best_test_m      <- test_m
        }
      }

      if (!is.na(best_tag)) {
        cv_results_gc[[length(cv_results_gc) + 1]] <- data.frame(
          species       = sp,
          condition     = cond,
          cv_method     = "LOMO",
          fold          = paste0("test_month_", test_month),
          n_train       = length(train_dates),
          n_test        = length(test_dates),
          best_tag      = best_tag,
          train_r2      = best_train_m$r_squared,
          train_kge     = best_train_m$kge,
          train_rmse    = best_train_m$rmse,
          train_bias    = best_train_m$bias,
          train_slope   = best_train_m$slope,
          train_composite = best_train_score,
          test_r2       = best_test_m$r_squared,
          test_kge      = best_test_m$kge,
          test_rmse     = best_test_m$rmse,
          test_bias     = best_test_m$bias,
          test_slope    = best_test_m$slope,
          test_composite = compute_composite(
            best_test_m$r_squared, best_test_m$kge,
            best_test_m$slope, best_test_m$bias, best_test_m$nrmse_pct
          )
        )
      }
    }
  }
}

cv_gc_df <- bind_rows(cv_results_gc)

# ==========================================================================
# 7. CROSS-VALIDATION: TEMPORAL 80/20 BLOCK SPLIT for Gc
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  CROSS-VALIDATION 2: TEMPORAL 80/20 BLOCK SPLIT (Gc)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cv_block_gc <- list()

for (sp in names(parameter_bases)) {
  sp_low <- species_lower_map[sp]

  for (cond in c("control", "drought")) {

    obs_sub <- obs_gc %>%
      filter(species == sp_low, treatment == cond) %>%
      group_by(date) %>%
      summarise(G_obs = mean(G_ms, na.rm = TRUE), .groups = "drop") %>%
      filter(month(date) %in% growing_months) %>%
      arrange(date)

    if (nrow(obs_sub) < 20) next

    n_total <- nrow(obs_sub)
    split_idx <- floor(n_total * 0.8)

    # 80% train (early), 20% test (late)
    train_dates <- obs_sub$date[1:split_idx]
    test_dates  <- obs_sub$date[(split_idx + 1):n_total]

    # Also do reverse: train on late, test on early (robustness check)
    # This tests if seasonal patterns are consistent
    splits <- list(
      early_train = list(train = train_dates, test = test_dates, label = "train_early_test_late"),
      late_train  = list(train = test_dates, test = train_dates, label = "train_late_test_early")
    )

    for (spl in splits) {
      if (length(spl$train) < 10 || length(spl$test) < 5) next

      best_train_score <- -Inf
      best_tag         <- NA_character_

      for (tag in names(lpj_cache[[sp]][[cond]])) {
        lpj_data <- lpj_cache[[sp]][[cond]][[tag]]
        if (is.null(lpj_data)) next

        train_joined <- lpj_data %>%
          filter(date %in% spl$train) %>%
          inner_join(obs_sub %>% select(date, G_obs), by = "date") %>%
          filter(Gc <= 12)

        if (nrow(train_joined) < 5) next
        train_m <- compute_metrics(train_joined$Gc, train_joined$G_obs)
        train_score <- compute_composite(
          train_m$r_squared, train_m$kge,
          train_m$slope, train_m$bias, train_m$nrmse_pct
        )

        if (train_score > best_train_score) {
          best_train_score <- train_score
          best_tag         <- tag
          best_train_m     <- train_m

          test_joined <- lpj_data %>%
            filter(date %in% spl$test) %>%
            inner_join(obs_sub %>% select(date, G_obs), by = "date") %>%
            filter(Gc <= 12)
          best_test_m <- compute_metrics(test_joined$Gc, test_joined$G_obs)
        }
      }

      if (!is.na(best_tag)) {
        cv_block_gc[[length(cv_block_gc) + 1]] <- data.frame(
          species       = sp,
          condition     = cond,
          cv_method     = "block_80_20",
          fold          = spl$label,
          n_train       = length(spl$train),
          n_test        = length(spl$test),
          best_tag      = best_tag,
          train_r2      = best_train_m$r_squared,
          train_kge     = best_train_m$kge,
          train_rmse    = best_train_m$rmse,
          train_bias    = best_train_m$bias,
          train_slope   = best_train_m$slope,
          train_composite = best_train_score,
          test_r2       = best_test_m$r_squared,
          test_kge      = best_test_m$kge,
          test_rmse     = best_test_m$rmse,
          test_bias     = best_test_m$bias,
          test_slope    = best_test_m$slope,
          test_composite = compute_composite(
            best_test_m$r_squared, best_test_m$kge,
            best_test_m$slope, best_test_m$bias, best_test_m$nrmse_pct
          )
        )
      }
    }
  }
}

cv_block_gc_df <- bind_rows(cv_block_gc)

# ==========================================================================
# 8. CROSS-VALIDATION: LEAVE-ONE-TREATMENT-OUT (LOTO) for Gc
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  CROSS-VALIDATION 3: LEAVE-ONE-TREATMENT-OUT (Gc)\n")
cat("  Tests: params from control → predict drought? (and vice versa)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cv_loto_gc <- list()

for (sp in names(parameter_bases)) {
  sp_low <- species_lower_map[sp]

  for (train_cond in c("control", "drought")) {
    test_cond <- setdiff(c("control", "drought"), train_cond)

    train_obs <- obs_gc %>%
      filter(species == sp_low, treatment == train_cond) %>%
      group_by(date) %>%
      summarise(G_obs = mean(G_ms, na.rm = TRUE), .groups = "drop") %>%
      filter(month(date) %in% growing_months)

    test_obs <- obs_gc %>%
      filter(species == sp_low, treatment == test_cond) %>%
      group_by(date) %>%
      summarise(G_obs = mean(G_ms, na.rm = TRUE), .groups = "drop") %>%
      filter(month(date) %in% growing_months)

    if (nrow(train_obs) < 10 || nrow(test_obs) < 10) next

    best_train_score <- -Inf
    best_tag         <- NA_character_

    # Search over training condition's parameter space
    # NOTE: we use train_cond's parameter bases AND train_cond's LPJ output
    # because parameters are calibrated to the training condition
    for (tag in names(lpj_cache[[sp]][[train_cond]])) {
      lpj_data <- lpj_cache[[sp]][[train_cond]][[tag]]
      if (is.null(lpj_data)) next

      # BUT: evaluate on training condition's observations
      train_joined <- lpj_data %>%
        inner_join(train_obs %>% select(date, G_obs), by = "date") %>%
        filter(Gc <= 12)

      if (nrow(train_joined) < 5) next
      train_m <- compute_metrics(train_joined$Gc, train_joined$G_obs)
      train_score <- compute_composite(
        train_m$r_squared, train_m$kge,
        train_m$slope, train_m$bias, train_m$nrmse_pct
      )

      if (train_score > best_train_score) {
        best_train_score <- train_score
        best_tag         <- tag
        best_train_m     <- train_m

        # Now evaluate SAME parameter set but on test condition
        # Need to use train_cond's LPJ output (same parameters)
        # against test_cond's observations
        test_joined <- lpj_data %>%
          inner_join(test_obs %>% select(date, G_obs), by = "date") %>%
          filter(Gc <= 12)
        best_test_m <- compute_metrics(test_joined$Gc, test_joined$G_obs)
      }
    }

    if (!is.na(best_tag)) {
      cv_loto_gc[[length(cv_loto_gc) + 1]] <- data.frame(
        species       = sp,
        train_condition = train_cond,
        test_condition  = test_cond,
        cv_method     = "LOTO",
        fold          = paste0("train_", train_cond, "_test_", test_cond),
        n_train       = nrow(train_obs),
        n_test        = nrow(test_obs),
        best_tag      = best_tag,
        train_r2      = best_train_m$r_squared,
        train_kge     = best_train_m$kge,
        train_rmse    = best_train_m$rmse,
        train_bias    = best_train_m$bias,
        train_slope   = best_train_m$slope,
        train_composite = best_train_score,
        test_r2       = best_test_m$r_squared,
        test_kge      = best_test_m$kge,
        test_rmse     = best_test_m$rmse,
        test_bias     = best_test_m$bias,
        test_slope    = best_test_m$slope,
        test_composite = compute_composite(
          best_test_m$r_squared, best_test_m$kge,
          best_test_m$slope, best_test_m$bias, best_test_m$nrmse_pct
        )
      )
      cat(sprintf("  %-8s train=%s → test=%s: best=%s  train R²=%.3f  test R²=%.3f\n",
                  sp, train_cond, test_cond, best_tag,
                  best_train_m$r_squared, best_test_m$r_squared))
    }
  }
}

cv_loto_gc_df <- bind_rows(cv_loto_gc)

# ==========================================================================
# 9. CROSS-VALIDATION: LEAVE-ONE-OUT for PSI (sparse data)
#    With only ~6 obs per species×condition, we can only do LOO.
#    Also add pooled-species LOO for more robust assessment.
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  CROSS-VALIDATION 4: LEAVE-ONE-OUT (PSI, limited data)\n")
cat("  WARNING: Only ~6 PSI obs per species×treatment. Results have high variance.\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cv_loo_psi <- list()

for (sp in names(parameter_bases)) {
  for (cond in c("control", "drought")) {

    obs_sub <- obs_psi %>%
      filter(species == sp, treatment == cond) %>%
      arrange(date)

    if (nrow(obs_sub) < 5) {
      cat(sprintf("  [SKIP] %-8s %-10s: only %d PSI obs\n", sp, cond, nrow(obs_sub)))
      next
    }

    # Leave-one-out: each observation date is a test fold
    for (i in seq_len(nrow(obs_sub))) {
      test_date <- obs_sub$date[i]
      train_obs <- obs_sub[-i, ]

      best_train_score <- -Inf
      best_tag         <- NA_character_

      for (tag in names(lpj_cache[[sp]][[cond]])) {
        lpj_data <- lpj_cache[[sp]][[cond]][[tag]]
        if (is.null(lpj_data)) next

        train_joined <- lpj_data %>%
          inner_join(train_obs %>% select(date, psi_leaf_md_obs),
                     by = "date") %>%
          filter(!is.na(psi_leaf), !is.na(psi_leaf_md_obs))

        if (nrow(train_joined) < 3) next
        train_m <- compute_metrics(train_joined$psi_leaf,
                                    train_joined$psi_leaf_md_obs)
        train_score <- compute_composite(
          train_m$r_squared, train_m$kge,
          train_m$slope, train_m$bias, train_m$nrmse_pct
        )

        if (train_score > best_train_score) {
          best_train_score <- train_score
          best_tag         <- tag
          best_train_m     <- train_m

          test_row <- lpj_data %>%
            filter(date == test_date) %>%
            inner_join(obs_sub %>% filter(date == test_date) %>%
                       select(date, psi_leaf_md_obs), by = "date")
          # For single test point, compute simple error
          if (nrow(test_row) > 0) {
            best_test_err <- abs(test_row$psi_leaf[1] - test_row$psi_leaf_md_obs[1])
            best_test_psi_sim <- test_row$psi_leaf[1]
            best_test_psi_obs <- test_row$psi_leaf_md_obs[1]
          } else {
            best_test_err <- NA_real_
            best_test_psi_sim <- NA_real_
            best_test_psi_obs <- NA_real_
          }
        }
      }

      if (!is.na(best_tag)) {
        cv_loo_psi[[length(cv_loo_psi) + 1]] <- data.frame(
          species       = sp,
          condition     = cond,
          cv_method     = "LOO",
          fold          = paste0("test_", as.character(test_date)),
          n_train       = nrow(train_obs),
          n_test        = 1,
          best_tag      = best_tag,
          train_r2      = best_train_m$r_squared,
          train_kge     = best_train_m$kge,
          train_rmse    = best_train_m$rmse,
          train_bias    = best_train_m$bias,
          train_slope   = best_train_m$slope,
          train_composite = best_train_score,
          test_abs_error = best_test_err,
          test_psi_sim  = best_test_psi_sim,
          test_psi_obs  = best_test_psi_obs,
          test_date     = test_date
        )
      }
    }
  }
}

cv_loo_psi_df <- bind_rows(cv_loo_psi)

# ==========================================================================
# 10. COMPILE ALL RESULTS
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  CROSS-VALIDATION RESULTS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ---- Gc: LOMO Summary ----
if (nrow(cv_gc_df) > 0) {
  cat("--- Gc: Leave-One-Month-Out Summary ---\n")
  lomo_summary <- cv_gc_df %>%
    group_by(species, condition) %>%
    summarise(
      n_folds        = n(),
      mean_train_r2  = mean(train_r2, na.rm = TRUE),
      mean_train_kge = mean(train_kge, na.rm = TRUE),
      mean_test_r2   = mean(test_r2, na.rm = TRUE),
      mean_test_kge  = mean(test_kge, na.rm = TRUE),
      r2_drop        = mean_train_r2 - mean_test_r2,
      kge_drop       = mean_train_kge - mean_test_kge,
      .groups = "drop"
    )
  print(as.data.frame(lomo_summary))

  cat("\nInterpretation:\n")
  cat("  r2_drop < 0.1  → parameters generalize well across months\n")
  cat("  r2_drop > 0.3  → overfitting; parameters are month-specific\n")
  cat("  negative kge on test → model fails to predict held-out month\n\n")
}

# ---- Gc: Block 80/20 Summary ----
if (nrow(cv_block_gc_df) > 0) {
  cat("--- Gc: Temporal 80/20 Block Split Summary ---\n")
  block_summary <- cv_block_gc_df %>%
    group_by(species, condition, fold) %>%
    summarise(
      train_r2 = mean(train_r2, na.rm = TRUE),
      test_r2  = mean(test_r2, na.rm = TRUE),
      test_kge = mean(test_kge, na.rm = TRUE),
      .groups = "drop"
    )
  print(as.data.frame(block_summary))
  cat("\n")
}

# ---- Gc: LOTO Summary ----
if (nrow(cv_loto_gc_df) > 0) {
  cat("--- Gc: Leave-One-Treatment-Out Summary ---\n")
  print(as.data.frame(cv_loto_gc_df %>%
    select(species, train_condition, test_condition,
           best_tag, train_r2, test_r2, test_kge)))

  cat("\nInterpretation:\n")
  cat("  test_r2 close to train_r2 → parameters transfer across water regimes\n")
  cat("  test_r2 << train_r2 → parameters are regime-specific; model structure may need revision\n\n")
}

# ---- PSI: LOO Summary ----
if (nrow(cv_loo_psi_df) > 0) {
  cat("--- PSI: Leave-One-Out Summary (high variance — interpret cautiously) ---\n")
  psi_loo_summary <- cv_loo_psi_df %>%
    group_by(species, condition) %>%
    summarise(
      n_folds          = n(),
      mean_train_r2    = mean(train_r2, na.rm = TRUE),
      mean_abs_error   = mean(test_abs_error, na.rm = TRUE),
      sd_abs_error     = sd(test_abs_error, na.rm = TRUE),
      .groups = "drop"
    )
  print(as.data.frame(psi_loo_summary))

  cat("\nNote: With only ~6 PSI observations per species×treatment, these LOO\n")
  cat("results have HIGH VARIANCE. Each fold trains on only 4-5 points.\n")
  cat("Treat these as indicative, not conclusive. Consider pooling across\n")
  cat("species or using multi-year PSI data (2018-2025 available) for more\n")
  cat("robust validation.\n\n")
}

# ==========================================================================
# 11. PARAMETER SELECTION STABILITY
#     Do different CV folds select the SAME best parameters?
# ==========================================================================
cat("--- Parameter Selection Stability ---\n")
cat("(How often is the same parameter selected as 'best' across folds?)\n\n")

all_gc_cv <- bind_rows(
  cv_gc_df %>% mutate(cv_source = "LOMO"),
  cv_block_gc_df %>% mutate(cv_source = "block_80_20")
)

if (nrow(all_gc_cv) > 0) {
  stability <- all_gc_cv %>%
    group_by(species, condition) %>%
    summarise(
      n_folds             = n(),
      unique_best_tags    = n_distinct(best_tag),
      top_tag             = names(sort(table(best_tag), decreasing = TRUE))[1],
      top_tag_frequency   = sort(table(best_tag), decreasing = TRUE)[1] / n(),
      .groups = "drop"
    ) %>%
    mutate(
      stability_assessment = case_when(
        top_tag_frequency >= 0.7 ~ "HIGH (same params consistently best)",
        top_tag_frequency >= 0.4 ~ "MODERATE (some agreement)",
        top_tag_frequency < 0.4  ~ "LOW (params depend heavily on fold)"
      )
    )
  print(as.data.frame(stability))
  cat("\n")
}

# ==========================================================================
# 12. VISUALIZATION
# ==========================================================================

# --- 12a. Train vs Test R² scatter (Gc LOMO) ---
if (nrow(cv_gc_df) > 0) {
  p <- ggplot(cv_gc_df, aes(x = train_r2, y = test_r2, color = species, shape = condition)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      title = "Cross-Validation: Train vs Test R² (Gc, Leave-One-Month-Out)",
      subtitle = "Points near the 1:1 line → parameters generalize well",
      x = "R² (training months)", y = "R² (held-out test month)"
    ) +
    scale_color_manual(values = c("Oak"="#E69F00", "Beech"="#0072B2",
                                   "Spruce"="#009E73", "Pine"="#F0E442")) +
    xlim(0, 1) + ylim(0, 1) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(OUTPUT_DIR, "cv_Gc_train_vs_test_R2.png"),
         p, width = 9, height = 7, dpi = 300)

  # --- 12b. R² drop by species × condition ---
  p2 <- ggplot(cv_gc_df, aes(x = species, y = train_r2 - test_r2, fill = species)) +
    geom_boxplot() +
    facet_wrap(~ condition) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 0.3, linetype = "dotted", color = "red") +
    labs(
      title = "R² Drop from Training to Test (Gc, LOMO)",
      subtitle = "Small drop → robust. Large drop → overfitting.",
      y = "R²(train) - R²(test)", x = ""
    ) +
    scale_fill_manual(values = c("Oak"="#E69F00", "Beech"="#0072B2",
                                  "Spruce"="#009E73", "Pine"="#F0E442")) +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(file.path(OUTPUT_DIR, "cv_Gc_R2_drop.png"),
         p2, width = 8, height = 5, dpi = 300)
}

# --- 12c. LOTO: cross-treatment transfer ---
if (nrow(cv_loto_gc_df) > 0) {
  p3 <- ggplot(cv_loto_gc_df,
               aes(x = train_condition, y = test_r2, fill = test_condition)) +
    geom_col(position = "dodge") +
    facet_wrap(~ species) +
    labs(
      title = "Cross-Treatment Transfer: Test R² (Gc, LOTO)",
      subtitle = "Parameters calibrated on one treatment, tested on the other",
      x = "Training condition", y = "Test R²"
    ) +
    scale_fill_manual(values = c("control" = "steelblue", "drought" = "coral")) +
    theme_minimal()

  ggsave(file.path(OUTPUT_DIR, "cv_Gc_cross_treatment.png"),
         p3, width = 9, height = 7, dpi = 300)
}

# --- 12d. Parameter Stability Heatmap ---
if (nrow(all_gc_cv) > 0) {
  # How often is each parameter tag selected as best?
  tag_freq <- all_gc_cv %>%
    group_by(species, condition) %>%
    mutate(total_folds = n()) %>%
    group_by(species, condition, best_tag, total_folds) %>%
    summarise(frequency = n() / first(total_folds), .groups = "drop") %>%
    group_by(species, condition) %>%
    slice_max(frequency, n = 5)  # top 5 per species×condition

  p4 <- ggplot(tag_freq, aes(x = reorder(best_tag, frequency), y = frequency, fill = species)) +
    geom_col() +
    facet_grid(species ~ condition, scales = "free") +
    coord_flip() +
    labs(
      title = "Parameter Selection Frequency Across CV Folds",
      subtitle = "Higher = same parameters chosen consistently → more robust",
      x = "Parameter combination", y = "Selection frequency"
    ) +
    scale_fill_manual(values = c("Oak"="#E69F00", "Beech"="#0072B2",
                                  "Spruce"="#009E73", "Pine"="#F0E442")) +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(file.path(OUTPUT_DIR, "cv_parameter_stability.png"),
         p4, width = 10, height = 8, dpi = 300)
}

# ==========================================================================
# 13. EXPORT CV RESULTS
# ==========================================================================
if (nrow(cv_gc_df) > 0) {
  write.csv(cv_gc_df, file.path(OUTPUT_DIR, "cv_Gc_LOMO.csv"), row.names = FALSE)
}
if (nrow(cv_block_gc_df) > 0) {
  write.csv(cv_block_gc_df, file.path(OUTPUT_DIR, "cv_Gc_block_80_20.csv"), row.names = FALSE)
}
if (nrow(cv_loto_gc_df) > 0) {
  write.csv(cv_loto_gc_df, file.path(OUTPUT_DIR, "cv_Gc_LOTO.csv"), row.names = FALSE)
}
if (nrow(cv_loo_psi_df) > 0) {
  write.csv(cv_loo_psi_df, file.path(OUTPUT_DIR, "cv_PSI_LOO.csv"), row.names = FALSE)
}

# Combined Gc CV summary
all_gc_cv_export <- bind_rows(
  cv_gc_df %>% mutate(cv_source = "LOMO"),
  cv_block_gc_df %>% mutate(cv_source = "block_80_20"),
  cv_loto_gc_df %>% mutate(cv_source = "LOTO")
)
if (nrow(all_gc_cv_export) > 0) {
  write.csv(all_gc_cv_export, file.path(OUTPUT_DIR, "cv_Gc_all_methods.csv"), row.names = FALSE)
}

# ==========================================================================
# 14. FINAL VERDICT
# ==========================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  CROSS-VALIDATION VERDICT\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (nrow(cv_gc_df) > 0) {
  overall_r2_drop <- mean(cv_gc_df$train_r2 - cv_gc_df$test_r2, na.rm = TRUE)
  overall_test_kge <- mean(cv_gc_df$test_kge, na.rm = TRUE)

  cat(sprintf("  Gc LOMO: mean test KGE = %.3f, mean R² drop = %.3f\n",
              overall_test_kge, overall_r2_drop))

  if (overall_r2_drop < 0.1 && overall_test_kge > 0.5) {
    cat("  ✅ Model is ROBUST for Gc — parameters generalize well across months.\n")
  } else if (overall_r2_drop < 0.2 && overall_test_kge > 0.3) {
    cat("  ⚠️  Model is MODERATELY robust for Gc — some overfitting, but reasonable.\n")
  } else {
    cat("  ❌ Model shows SIGNIFICANT OVERFITTING for Gc — parameters do NOT generalize.\n")
    cat("     Consider: (1) reducing parameter space, (2) multi-site calibration,\n")
    cat("     (3) adding prior constraints on parameters.\n")
  }
}

if (nrow(cv_loto_gc_df) > 0) {
  loto_test_kge <- mean(cv_loto_gc_df$test_kge, na.rm = TRUE)
  cat(sprintf("\n  Gc LOTO: mean cross-treatment test KGE = %.3f\n", loto_test_kge))
  if (loto_test_kge > 0.3) {
    cat("  ✅ Parameters transfer reasonably between control and drought.\n")
  } else {
    cat("  ⚠️  Parameters do NOT transfer well between treatments.\n")
    cat("     The hydraulic parameterization may need to differ by water regime.\n")
  }
}

if (nrow(cv_loo_psi_df) > 0) {
  cat(sprintf("\n  PSI LOO: mean absolute error = %.3f MPa (n=%d folds)\n",
              mean(cv_loo_psi_df$test_abs_error, na.rm = TRUE),
              nrow(cv_loo_psi_df)))
  cat("  ⚠️  PSI validation is limited by sparse observations (6 per species×treatment).\n")
  cat("  Recommend: include multi-year PSI data (2018-2025) for more robust validation.\n")
}

cat("\n*** Cross-validation complete. Results saved to:", OUTPUT_DIR, "***\n")