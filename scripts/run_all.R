# Entry point to run the full analysis end-to-end.
#
# Run from the REPOSITORY ROOT, not from scripts/:
#   Rscript scripts/run_all.R
# All path constants in scripts/utils_functions.R (DIR_INPUT <- "data",
# results/tables, results/figures, results/rds, results/stan) are relative
# to the repository root.
source("scripts/utils_functions.R"); ensure_packages()

# Runtime: steps 01-04, 06, 08, 11 and 15 are quick (06 only defines the
# growth-model functions). Step 07 can take a while but does not use Stan: it fits five
# candidate growth models to every monoculture time series by multi-start
# nonlinear least squares, and regenerates the per-taxon curve-fit plots and
# parameter tables behind the supplementary growth-curve figures. Note that 07
# clears the workspace (rm(list = ls())); every later step re-sources its own
# helpers, so this is harmless. Steps that fit Stan models are expensive,
# although each caches its fits in results/rds/ and reloads them on rerun
# (delete a cached .rds to force that refit):
#   05 (main fit and blocked cross-validation), 09 (leave-one-stress-out, 8 fits),
#   10 (alternative traits and log-linear growth rates), 12 (ablation, 15 fits),
#   13 (leave-one-taxon-out, 12 fits), 14 (singles to combinations, 1 fit).
# 08 reuses the cached main fit from 05 and does not refit.
steps <- c(
  "scripts/01_data_preprocessing.R",
  "scripts/02_fig1_diversity_permanova.R",
  "scripts/03_fig1_indicator_species.R",
  "scripts/04_fig2_growth_zscores.R",
  "scripts/05_fig3_bayes_od.R",
  "scripts/06_growth_curve_models.R",
  "scripts/07_growth_curves.R",
  "scripts/08_biomass_decomposition.R",
  "scripts/09_loso_cv.R",
  "scripts/10_trait_sensitivity.R",
  "scripts/11_winner_takes_all.R",
  "scripts/12_composition_ablation.R",
  "scripts/13_leave_one_taxon_out.R",
  "scripts/14_singles_to_combos.R",
  "scripts/15_main_figures.R"
)
for (s in steps) {
  message(">>> Running: ", s)
  source(s, local = TRUE)
}

# Save session info
try({
  sink(file = file.path("results", "sessionInfo.txt"))
  print(Sys.time()); print(R.version.string); sessionInfo()
  sink()
}, silent = TRUE)
message("All steps completed. Outputs in results/")
