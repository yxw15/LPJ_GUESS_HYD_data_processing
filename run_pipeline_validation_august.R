# ==========================================================================
# AUGUST VALIDATION PIPELINE
# Runs all plot scripts filtered to August ONLY.
# Used to validate parameters trained on July+September.
#
# Usage (after running LPJ-GUESS with best July+Sept-trained parameters):
#   Rscript run_pipeline_validation_august.R
#
# Before running, update the LPJ output path if it differs from the default.
# ==========================================================================

script_folder <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/R_scripts"

# August-only pipeline: same scripts, but each filters to month 8
scripts_to_run <- c(
  "tidy_lpj_twd_output_daily_august.R",
  "plot_compare_Gc_lpjtwd_hoelstein_august.R",
  "plot_compare_PSI_lpjtwd_hoelstein_august.R",
  "plot_compare_twd_lpjtwd_hoelstein_august.R",
  "plot_lpj_stem_storage_Gc_ET_TWD_vpd_august.R",
  "plot_compare_Gc_PsiL_lpjtwd_hoelstein_august.R",
  "plot_compare_wcont_lpjtwd_hoelstein_august.R"
)

for (file in scripts_to_run) {
  full_path <- file.path(script_folder, file)
  if (file.exists(full_path)) {
    message("Running: ", file)
    source(full_path, local = new.env())
  } else {
    message("SKIPPING (file not found): ", file)
  }
}

message("\n=== August validation pipeline complete ===")
message("Figures saved to Figures/lpj_guess_stem_storage/validation_august/")
