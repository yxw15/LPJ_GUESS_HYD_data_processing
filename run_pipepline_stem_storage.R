script_folder <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/R_scripts"

# List your files here
scripts_to_run <- c(
  "tidy_lpj_twd_output_daily.R",
  "plot_compare_Gc_lpjtwd_hoelstein.R",
  "plot_compare_PSI_lpjtwd_hoelstein.R",
  "plot_compare_twd_lpjtwd_hoelstein.R",
  "plot_lpj_stem_storage_Gc_ET_TWD_vpd.R",
  "plot_compare_Gc_PsiL_lpjtwd_hoelstein.R",
  "plot_compare_twd_lpjtwd_hoelstein.R",
  "plot_compare_wcont_lpjtwd_hoelstein.R"
)

# Loop through and source each one
for (file in scripts_to_run) {
  # This creates the correct path safely
  full_path <- file.path(script_folder, file)
  
  message("Running: ", file)
  source(full_path, local = new.env())
}