# ============================================================
# RUN LPJ-GUESS
# Native scale = raw
# Coarser scale from daily/monthly = mean/min/max
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
})


# ============================================================
# Build daily raw + aggregated
# ============================================================

daily_all <- build_daily_all(
  base_dir = base_dir,
  daily_vars = daily_vars,
  species_map = species_map,
  year_min = year_min,
  year_max = year_max,
  day_starts_at_zero = day_starts_at_zero
)

daily_agg <- make_daily_aggregates(
  daily_df = daily_all,
  species_map = species_map,
  stats = stats_to_make
)

# ============================================================
# Daily-origin plots
# raw daily + aggregated monthly/yearly
# ============================================================

daily_plots <- list()

for (i in seq_len(nrow(daily_vars))) {
  stem_i      <- daily_vars$stem[i]
  var_label_i <- daily_vars$var_label[i]
  ylab_i      <- daily_vars$ylab[i]
  
  daily_plots[[stem_i]] <- list()
  
  # raw daily, ALL
  daily_plots[[stem_i]][["daily_raw_ALL"]] <- plot_raw_ts(
    df = daily_all,
    cb_palette = cb_palette,
    out_dir = out_dir,
    var_label = var_label_i,
    ylab = ylab_i,
    filename_stub = stem_i,
    scale_name = "daily",
    species = NULL,
    plot_year_min = plot_year_min,
    plot_year_max = plot_year_max
  )
  
  # monthly/yearly aggregated, ALL
  for (st in stats_to_make) {
    daily_plots[[stem_i]][[paste0("monthly_", st, "_ALL")]] <- plot_aggregated_ts(
      df = daily_agg,
      cb_palette = cb_palette,
      out_dir = out_dir,
      var_label = var_label_i,
      ylab = ylab_i,
      filename_stub = stem_i,
      scale = "monthly",
      stat = st,
      species = NULL,
      plot_year_min = plot_year_min,
      plot_year_max = plot_year_max
    )
    
    daily_plots[[stem_i]][[paste0("yearly_", st, "_ALL")]] <- plot_aggregated_ts(
      df = daily_agg,
      cb_palette = cb_palette,
      out_dir = out_dir,
      var_label = var_label_i,
      ylab = ylab_i,
      filename_stub = stem_i,
      scale = "yearly",
      stat = st,
      species = NULL,
      plot_year_min = plot_year_min,
      plot_year_max = plot_year_max
    )
  }
  
  for (sp in species_map$species) {
    daily_plots[[stem_i]][[paste0("daily_raw_", sp)]] <- plot_raw_ts(
      df = daily_all,
      cb_palette = cb_palette,
      out_dir = out_dir,
      var_label = var_label_i,
      ylab = ylab_i,
      filename_stub = stem_i,
      scale_name = "daily",
      species = sp,
      plot_year_min = plot_year_min,
      plot_year_max = plot_year_max
    )
    
    for (st in stats_to_make) {
      daily_plots[[stem_i]][[paste0("monthly_", st, "_", sp)]] <- plot_aggregated_ts(
        df = daily_agg,
        cb_palette = cb_palette,
        out_dir = out_dir,
        var_label = var_label_i,
        ylab = ylab_i,
        filename_stub = stem_i,
        scale = "monthly",
        stat = st,
        species = sp,
        plot_year_min = plot_year_min,
        plot_year_max = plot_year_max
      )
      
      daily_plots[[stem_i]][[paste0("yearly_", st, "_", sp)]] <- plot_aggregated_ts(
        df = daily_agg,
        cb_palette = cb_palette,
        out_dir = out_dir,
        var_label = var_label_i,
        ylab = ylab_i,
        filename_stub = stem_i,
        scale = "yearly",
        stat = st,
        species = sp,
        plot_year_min = plot_year_min,
        plot_year_max = plot_year_max
      )
    }
  }
}

# ============================================================
# Monthly-origin plots
# raw monthly + aggregated yearly
# ============================================================

monthly_plots <- list()

if (nrow(monthly_vars) > 0) {
  for (i in seq_len(nrow(monthly_vars))) {
    file_i      <- monthly_vars$file[i]
    var_label_i <- monthly_vars$var_label[i]
    ylab_i      <- monthly_vars$ylab[i]
    stem_i      <- monthly_vars$stem[i]
    
    monthly_raw <- build_monthly_all(
      base_dir = base_dir,
      filename = file_i,
      var_label = var_label_i,
      species_map = species_map,
      year_min = year_min,
      year_max = year_max
    )
    
    monthly_agg <- make_monthly_aggregates(
      monthly_df = monthly_raw,
      species_map = species_map,
      stats = stats_to_make
    )
    
    monthly_plots[[stem_i]] <- list()
    
    monthly_plots[[stem_i]][["monthly_raw_ALL"]] <- plot_raw_ts(
      df = monthly_raw,
      cb_palette = cb_palette,
      out_dir = out_dir,
      var_label = var_label_i,
      ylab = ylab_i,
      filename_stub = stem_i,
      scale_name = "monthly",
      species = NULL,
      plot_year_min = plot_year_min,
      plot_year_max = plot_year_max
    )
    
    for (st in stats_to_make) {
      monthly_plots[[stem_i]][[paste0("yearly_", st, "_ALL")]] <- plot_aggregated_ts(
        df = monthly_agg,
        cb_palette = cb_palette,
        out_dir = out_dir,
        var_label = var_label_i,
        ylab = ylab_i,
        filename_stub = stem_i,
        scale = "yearly",
        stat = st,
        species = NULL,
        plot_year_min = plot_year_min,
        plot_year_max = plot_year_max
      )
    }
    
    for (sp in species_map$species) {
      monthly_plots[[stem_i]][[paste0("monthly_raw_", sp)]] <- plot_raw_ts(
        df = monthly_raw,
        cb_palette = cb_palette,
        out_dir = out_dir,
        var_label = var_label_i,
        ylab = ylab_i,
        filename_stub = stem_i,
        scale_name = "monthly",
        species = sp,
        plot_year_min = plot_year_min,
        plot_year_max = plot_year_max
      )
      
      for (st in stats_to_make) {
        monthly_plots[[stem_i]][[paste0("yearly_", st, "_", sp)]] <- plot_aggregated_ts(
          df = monthly_agg,
          cb_palette = cb_palette,
          out_dir = out_dir,
          var_label = var_label_i,
          ylab = ylab_i,
          filename_stub = stem_i,
          scale = "yearly",
          stat = st,
          species = sp,
          plot_year_min = plot_year_min,
          plot_year_max = plot_year_max
        )
      }
    }
  }
}

# ============================================================
# Yearly-origin plots
# raw yearly only
# ============================================================

yearly_plots <- list()

if (nrow(yearly_vars) > 0) {
  for (i in seq_len(nrow(yearly_vars))) {
    file_i      <- yearly_vars$file[i]
    var_label_i <- yearly_vars$var_label[i]
    ylab_i      <- yearly_vars$ylab[i]
    stem_i      <- yearly_vars$stem[i]
    
    yearly_raw <- build_yearly_all(
      base_dir = base_dir,
      filename = file_i,
      var_label = var_label_i,
      species_map = species_map,
      year_min = year_min,
      year_max = year_max
    )
    
    yearly_plots[[stem_i]] <- list()
    
    yearly_plots[[stem_i]][["yearly_raw_ALL"]] <- plot_raw_ts(
      df = yearly_raw,
      cb_palette = cb_palette,
      out_dir = out_dir,
      var_label = var_label_i,
      ylab = ylab_i,
      filename_stub = stem_i,
      scale_name = "yearly",
      species = NULL,
      plot_year_min = plot_year_min,
      plot_year_max = plot_year_max
    )
    
    for (sp in species_map$species) {
      yearly_plots[[stem_i]][[paste0("yearly_raw_", sp)]] <- plot_raw_ts(
        df = yearly_raw,
        cb_palette = cb_palette,
        out_dir = out_dir,
        var_label = var_label_i,
        ylab = ylab_i,
        filename_stub = stem_i,
        scale_name = "yearly",
        species = sp,
        plot_year_min = plot_year_min,
        plot_year_max = plot_year_max
      )
    }
  }
}