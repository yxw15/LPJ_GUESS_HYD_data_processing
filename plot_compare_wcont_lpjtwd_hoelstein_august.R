VALIDATION_MONTH <- 8  # August only
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(ggplot2)

# set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# define and create the output directory for figures if it doesn't exist
output_dir = "Figures/lpj_guess_stem_storage/validation_august/Wcont"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==========================================
# 1. read and process model data (mwcont) from 10cm, 40cm, 80cm layers
# ==========================================
file_list = list.files(
  path = "results_lpj/results_hoelstein_stem_storage",
  pattern = "mwcont_(10cm|40cm|80cm)\\.out$",
  recursive = TRUE,
  full.names = TRUE
)

cat("number of mwcont files found:", length(file_list), "\n")
cat("files:\n")
for (f in file_list) cat("  ", f, "\n")

mwcont_all_combined = file_list %>%
  map_df(function(file_path) {
    path_parts = unlist(strsplit(file_path, "/"))

    # path_parts maps out to:
    # [1] "results_lpj" [2] "results_hoelstein_stem_storage" [3] treatment [4] species [5] file
    current_treatment = tolower(path_parts[3])
    current_species   = tolower(path_parts[4])
    # extract layer depth from filename (e.g. "10cm", "40cm", "80cm")
    current_layer     = str_extract(file_path, "10cm|40cm|80cm")

    read.table(file_path, header = TRUE) %>%
      pivot_longer(cols = Jan:Dec, names_to = "Month", values_to = "wcont_value") %>%
      mutate(Month = tolower(Month)) %>%
      unite("year-month", Year, Month, sep = "-") %>%
      mutate(
        wcont_value = wcont_value * 100,
        species     = current_species,
        treatment   = current_treatment,
        soil_layer  = current_layer
      ) %>%
      select(`year-month`, wcont_value, species, treatment, soil_layer)
  })

# validate data was loaded correctly
stopifnot(nrow(mwcont_all_combined) > 0)
stopifnot("soil_layer" %in% colnames(mwcont_all_combined))

cat("\nunique treatments found:", unique(mwcont_all_combined$treatment), "\n")
cat("unique species found:", unique(mwcont_all_combined$species), "\n")
cat("unique soil layers found:", unique(mwcont_all_combined$soil_layer), "\n")

# ==========================================
# 2. read and aggregate field data (sm) to monthly
# ==========================================
sm_raw = readRDS("SCCII/Soil_water_content_and_plant_water_uptake_depth/SM_REW_SEL_2018_2025.rds")
sm_raw = sm_raw %>%
  mutate(treatment = recode(treatment, "treatment" = "drought"))

sm_monthly = sm_raw %>%
  drop_na() %>%
  mutate(
    `year-month` = paste(year(date), tolower(as.character(month(date, label = TRUE))), sep = "-")
  ) %>%
  group_by(`year-month`, depth, treatment) %>%
  summarize(sm_mean = mean(SM, na.rm = TRUE), .groups = "drop")

# ==========================================
# 3. intersect & format date for all graphing
# ==========================================
common_months = intersect(mwcont_all_combined$`year-month`, sm_monthly$`year-month`)

cat("number of overlapping months found:", length(common_months), "\n")

mwcont_filtered = mwcont_all_combined %>% filter(`year-month` %in% common_months)
sm_filtered     = sm_monthly     %>% filter(`year-month` %in% common_months)

to_date_format = function(df) {
  df %>% mutate(plot_date = ym(str_replace(`year-month`, "-", " ")))
}

mwcont_filtered = to_date_format(mwcont_filtered)
sm_filtered     = to_date_format(sm_filtered)

# ==========================================
# 4. plots 1 & 2 (individual time series)
# ==========================================

# plot 1: modelled water content (mwcont) — faceted by species x treatment, colored by soil layer
p1 = ggplot(mwcont_filtered, aes(x = plot_date, y = wcont_value, color = soil_layer)) +
  geom_line(linewidth = 0.8) +
  facet_grid(species ~ treatment) +
  labs(
    title = "modelled soil water content (mwcont) over time",
    x = "date",
    y = "water content (%)",
    color = "soil layer"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%b", date_breaks = "6 months") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave(filename = file.path(output_dir, "plot1_modelled_mwcont.png"), plot = p1, width = 10, height = 7, dpi = 300)

# plot 2: measured soil moisture (sm)
p2 = ggplot(sm_filtered, aes(x = plot_date, y = sm_mean, color = as.factor(depth))) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ treatment) +
  labs(
    title = "field measured soil moisture (sm) monthly averages",
    x = "date",
    y = "sm mean (%)",
    color = "depth (cm)"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%b", date_breaks = "6 months") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
ggsave(filename = file.path(output_dir, "plot2_measured_sm.png"), plot = p2, width = 10, height = 5, dpi = 300)

# ==========================================
# 5. plot 3: combined model vs measured data (2 rows of treatments, grouped by species columns)
# ==========================================

# 1. Pull out the exact timeframe for the observed field drought data
drought_time_range = sm_filtered %>%
  filter(treatment == "drought") %>%
  summarize(
    min_date = min(plot_date, na.rm = TRUE),
    max_date = max(plot_date, na.rm = TRUE)
  )

min_drought_date = drought_time_range$min_date
max_drought_date = drought_time_range$max_date

# 2. prepare modelled data (filtered by drought window boundaries)
mw_to_merge = mwcont_filtered %>%
  filter(plot_date >= min_drought_date & plot_date <= max_drought_date) %>%
  mutate(
    value = wcont_value,
    data_source = "modelled (mwcont)",
    source_layer = case_when(
      soil_layer == "10cm" ~ "model: 10 cm",
      soil_layer == "40cm" ~ "model: 40 cm",
      soil_layer == "80cm" ~ "model: 80 cm"
    )
  ) %>%
  select(plot_date, value, species, treatment, data_source, source_layer)

# 3. prepare observed field data (filtered by drought window boundaries)
#    note: observed soil moisture is per-plot (not per-species), so we replicate
#    the observed data for each modelled species to enable direct comparison
unique_species = unique(mwcont_filtered$species)

sm_to_merge = unique_species %>%
  map_df(function(sp) {
    sm_filtered %>%
      filter(plot_date >= min_drought_date & plot_date <= max_drought_date) %>%
      mutate(
        species = sp,
        value = sm_mean,
        data_source = "measured (sm)",
        source_layer = case_when(
          depth == "10" ~ "observed: 10 cm",
          depth == "40" ~ "observed: 40 cm",
          depth == "80" ~ "observed: 80 cm",
          TRUE          ~ paste0("observed: ", depth, " cm")
        )
      )
  }) %>%
  select(plot_date, value, species, treatment, data_source, source_layer)

# 4. combine both frames
combined_data = bind_rows(mw_to_merge, sm_to_merge)

# Ensure species names are capitalized cleanly for the plot facets
combined_data = combined_data %>%
  mutate(species = str_to_title(species))

# 5. define custom distinct color hex codes (3 model layers + 3 observed depths)
custom_colors = c(
  "model: 10 cm"    = "#e66101",
  "model: 40 cm"    = "#fdb863",
  "model: 80 cm"    = "#b2182b",
  "observed: 10 cm" = "#4d9221",
  "observed: 40 cm" = "#2b8cbe",
  "observed: 80 cm" = "#045a8d"
)

# 6. final plot execution (Rows = treatment, Columns = species)
p3 = ggplot(combined_data, aes(x = plot_date, y = value, color = source_layer, linetype = data_source)) +
  geom_line(linewidth = 0.9) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("modelled (mwcont)" = "solid", "measured (sm)" = "dashed")) +
  labs(
    title = "validation: modelled (mwcont) vs field observed (sm) soil water dynamics",
    subtitle = paste("synchronized to observed drought experimental window:", min_drought_date, "to", max_drought_date),
    x = "date",
    y = "volumetric water content / soil moisture (%)",
    color = "data streams & depths",
    linetype = "data source"
  ) +
  theme_minimal(base_size = 11) +
  scale_x_date(date_labels = "%Y-%b", date_breaks = "3 months") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p3)
ggsave(filename = file.path(output_dir, "plot3_combined_validation.png"), plot = p3, width = 15, height = 8, dpi = 300)

# ==========================================
# 6. plot 4: distribution of values (boxplot) by treatment and species
# ==========================================
p4 = ggplot(combined_data, aes(x = source_layer, y = value, fill = source_layer)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_grid(treatment ~ species) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "distribution of modelled vs observed soil moisture values",
    subtitle = paste("analysis period:", min_drought_date, "to", max_drought_date),
    x = "data streams & depths",
    y = "volumetric water content / soil moisture (%)",
    fill = "data streams & depths"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p4)
ggsave(filename = file.path(output_dir, "plot4_distribution_boxplot.png"), plot = p4, width = 15, height = 8, dpi = 300)

# ==========================================
# 7. plot 5: distribution of values for August only (summer)
# ==========================================
# Filter down data exclusively to June (6), July (7), August (8), and September (9)
summer_data = combined_data %>%
  filter(month(plot_date) == VALIDATION_MONTH)

p5 = ggplot(summer_data, aes(x = source_layer, y = value, fill = source_layer)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_grid(treatment ~ species) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "distribution of soil moisture values (August only)",
    subtitle = paste("August comparison within window:", min_drought_date, "to", max_drought_date),
    x = "data streams & depths",
    y = "volumetric water content / soil moisture (%)",
    fill = "data streams & depths"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p5)
ggsave(filename = file.path(output_dir, "plot5_summer_distribution_boxplot.png"), plot = p5, width = 15, height = 8, dpi = 300)
