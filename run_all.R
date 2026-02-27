
# Entry point to run the full analysis end-to-end
source("scripts/utils_functions.R"); ensure_packages()

steps <- c(
  "scripts/01_data_preprocessing.R",
  "scripts/02_fig1_diversity_permanova.R",
  "scripts/03_fig1_indicator_species.R",
  "scripts/04_fig2_growth_zscores.R",
  "scripts/05_fig3_bayes_od.R"
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
