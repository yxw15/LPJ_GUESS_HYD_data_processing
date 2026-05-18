setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

psiL_hoelstein <- read.csv("SCCII/psiL_hoelstein.csv")

sap_daily <- read.csv("SCCII/sap_daily.csv")
sap_daily_control <- sap_daily %>% filter(treatment == "control")

# 1. Get the list of unique control tree IDs from your sap data
control_tree_ids <- sap_daily_control %>%
  distinct(tree_id)

# 2. Filter psiL data to only include these tree IDs
psiL_control <- psiL_hoelstein %>%
  semi_join(control_tree_ids, by = "tree_id")

# 3. Save the result
write_csv(psiL_control, "SCCII/psiL_hoelstein_control.csv")
write_csv(sap_daily_control, "SCCII/sap_hoelstein_control.csv")

# 4. Optional: Check how many trees were kept
cat("Number of unique control trees found in psiL:", 
    length(unique(psiL_control$tree_id)), "\n")

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

psiL_hoelstein <- read.csv("SCCII/psiL_hoelstein.csv")

#---------------------------------------------------

psiL_hoelstein <- read.csv("SCCII/psiL_hoelstein.csv")
sap_daily <- read.csv("SCCII/sap_daily.csv")
sap_daily_drought <- sap_daily %>% filter(treatment == "treatment")

# 1. Get the list of unique control tree IDs from your sap data
drought_tree_ids <- sap_daily_drought %>%
  distinct(tree_id)

# 2. Filter psiL data to only include these tree IDs
psiL_drought <- psiL_hoelstein %>%
  semi_join(drought_tree_ids, by = "tree_id")

# 3. Save the result
write_csv(psiL_drought, "SCCII/psiL_hoelstein_drought.csv")
write_csv(sap_daily_drought, "SCCII/sap_hoelstein_drought.csv")

# 4. Optional: Check how many trees were kept
cat("Number of unique drought trees found in psiL:", 
    length(unique(psiL_drought$tree_id)), "\n")
