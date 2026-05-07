setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

psiS_hoelstein_raw <- readRDS("SCCII/soilwaterpotential/Hoelstein_SWP_2023_2025_ALL_rescaled.rds")

# 2. Filter for control and calculate daily means
psiS_hoelstein_daily <- psiS_hoelstein_raw %>%
  filter(treatment == "control") %>%
  mutate(date = as.Date(timestamp_UTC)) %>%
  group_by(date) %>%
  summarise(
    # Calculating the site-level mean for the rescaled SWP
    psiS_mean = mean(SWP_site_rescaled, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Save as CSV
write_csv(psiS_hoelstein_daily, "SCCII/psiS_hoelstein_control.csv")