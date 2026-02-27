# Scaling the impacts of multiple stressors from populations to ecosystems

This repository contains the R scripts used for analyses and figures in the manuscript:
**Carmichael H., Cakin I., Busi S.B., Read D., Yvon-Durocher G.**

This README has been updated to match the current script set and output file names (including revised supplementary table/figure numbering used in the scripts).

## Repository structure (current scripts)

### Main analysis scripts
- `01_data_preprocessing.R` — builds and filters the phyloseq object and exports cleaned RDS objects.
- `02_fig1_diversity_permanova.R` — Shannon diversity, db-RDA, PERMANOVA pairwise heatmap (Fig 1a–c + tables).
- `03_fig1_indicator_species.R` — indicator species analysis and Fig 1d + full ISA table.
- `04_fig2_growth_zscores.R` — growth-rate harmonisation and z-score summaries (Fig 2 + tables).
- `05_fig3_bayes_od.R` — Bayesian composition model + OD prediction analyses (Fig 3 and supplementary outputs), plus blocked cross-validation.

### Helper scripts
- `utils_functions.R` — shared package loading, path helpers, palettes, plot labels, and utility functions.
- `growth_curve_models.R` — growth model definitions (Baranyi, Gompertz, Buchanan, etc.) and Sharpe–Schoolfield helper functions.

### Additional growth-curve fitting script (upstream / optional)
- `GrowthCurves.R` — standalone growth-curve fitting pipeline (AICc model selection across multiple growth models) used to generate growth-fit outputs for downstream analyses.

## Expected folders and path conventions

The scripts use path helpers in `utils_functions.R` and expect / create the following directories:

- `data/` — input files (CSV/TSV/tree files)
- `results/figures/` — figure outputs (`.tiff`)
- `results/tables/` — table outputs (`.csv`)
- `results/rds/` — intermediate and model objects (`.rds`)

> **Important:** The scripts call helper files as `source("scripts/utils_functions.R")` and `source("scripts/growth_curve_models.R")`. If your files are currently in the project root, either move them into a `scripts/` folder or update the `source()` paths.

## R version and package requirements

### R
- R >= 4.3 recommended

### Core packages used across the current pipeline
`utils_functions.R` installs / loads most requirements automatically via `ensure_packages()`.

Main packages include:
- `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggpubr`, `scales`, `ggrepel`, `patchwork`
- `ape`, `phyloseq`, `vegan`, `pairwiseAdonis`, `ggtree`, `stringr`, `picante`, `phangorn`, `microeco`, `reshape2`
- `broom`, `car`, `effectsize`, `emmeans`, `readr`, `permute`, `indicspecies`, `ragg`
- `rstan`, `loo`, `philentropy`

Additional packages used by `GrowthCurves.R`:
- `nls.multstart`, `MuMIn`, `purrr`, `readr`

## Analysis pipeline (current)

### 1) Data preprocessing

**Script:** `01_data_preprocessing.R`

**Inputs (in `data/`):**
- `ASV.tax.tsv`
- `ASV.counts.tsv`
- `Sequencing_metadata2.csv`
- `ASV1-13_tree.txt`

**What it does:**
- Loads ASV taxonomy, count matrix, sample metadata, and phylogenetic tree.
- Builds a `phyloseq` object.
- Converts to relative abundance and filters taxa by mean abundance (>0.5%).
- Optionally merges two *Pseudomonas* ASVs (if present).
- Subsets to `Diversity == "Eight"` and non-evolved communities (`Evo_status != "Evo"`).
- Aligns abundance and metadata tables by common samples.

**Outputs (`results/rds/`):**
- `phyloseq_filtered.rds`
- `ASV_matrix_clean.rds`
- `metadata_clean.rds`
- `taxa_table.rds`

---

### 2) Fig 1a–c: Diversity, db-RDA, and PERMANOVA

**Script:** `02_fig1_diversity_permanova.R`

**Inputs:**
- `results/rds/ASV_matrix_clean.rds`
- `results/rds/metadata_clean.rds`

**What it does:**
- Recodes stress labels (e.g. `Salinity` -> `Sal`, `Temperature` -> `Temp`).
- Computes Bray–Curtis distances.
- Runs db-RDA using pH, salinity, and temperature gradients.
- Runs global and pairwise PERMANOVA on community composition.
- Computes Shannon diversity, ANOVA, effect sizes, and Tukey post-hoc contrasts.

**Figure outputs (`results/figures/`):**
- `Fig_1a_Shannon_diversity_boxplot_with_signif.tiff`
- `Fig_1b_dbRDA_gradients.tiff`
- `Fig_1c_PERMANOVA_pairwise_heatmap.tiff`

**Table outputs (`results/tables/`):**
- `Table_S5_Tukey_posthoc_effectsizes.csv`
- `Table_S6_PERMANOVA_pairwise.csv`
- `Shannon_ANOVA_effectsizes.csv`

> Note: Compared with older versions, the script currently writes **Tukey post-hoc** as **Table S5** and **pairwise PERMANOVA** as **Table S6**.

---

### 3) Fig 1d: Indicator species analysis

**Script:** `03_fig1_indicator_species.R`

**Inputs:**
- `results/rds/ASV_matrix_clean.rds`
- `results/rds/metadata_clean.rds`
- `results/rds/taxa_table.rds`

**What it does:**
- Runs indicator species analysis (`indicspecies::multipatt`, `IndVal.g`) using one-vs-rest contrasts for each stress condition.
- Harmonises genus labels (including disambiguating *Pseudomonas* taxa).
- Extracts the top indicator taxon per stress condition (`TopN = 1`) and plots indicator values.

**Figure output (`results/figures/`):**
- `Fig_1d_ISA_Top1_1vsRest_singlepanel.tiff`

**Table output (`results/tables/`):**
- `Table_S7_ISA_stats_full.csv`

---

### 4) Fig 2: Growth z-scores across stress combinations

**Script:** `04_fig2_growth_zscores.R`

**Inputs (in `data/`):**
- `LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv` (**required**)
- `LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv` (optional if available)
- `LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv` (optional if available)

**What it does:**
- Loads growth-fit outputs from pH / temperature / salinity experiments.
- Harmonises columns (`Id`, growth rate `g`, replicate, and condition levels).
- Normalises pH / temperature / salinity levels to the experimental settings.
- Derives a unified `Stress` factor (Control, single stressors, and combinations).
- Builds a combined long table across sources and filters invalid/negative growth rates.
- Computes within-stress z-scores and summarises mean ± SE by taxon and stress.
- Produces the compact Fig 2 dotplot.

**Figure output (`results/figures/`):**
- `Fig_2_DOTPLOT_z_by_stress_compact.tiff`

**Table outputs (`results/tables/`):**
- `Table_S8a_Growth_all_sources_stressors_long.csv`
- `Table_S8_Zscores_summary.csv`

**RDS outputs (`results/rds/`):**
- `Growth_all_sources_stressors_long.rds`
- `Zscores_objects.rds`

---

### 5) Fig 3 + supplementary Bayesian composition / OD analysis + blocked CV

**Script:** `05_fig3_bayes_od.R`

This is a combined script with two parts:
- **Part A:** Main Bayesian composition model + OD predictions (Fig 3 and supplementary analyses)
- **Part B:** Blocked cross-validation (out-of-sample evaluation)

#### Part A inputs (in `data/`)
Core inputs include:
- `All_Abundance_Data.csv`
- `CommunityOD_All_Jan23_edit_zeros.csv`
- single-taxon growth files loaded inside the script (stress-specific growth-fit tables)

#### Part A overview
- Fits a **Dirichlet–softmax Bayesian model** (via `rstan`) linking z-scored growth traits to observed community composition.
- Generates posterior predicted compositions and abundance-weighted mean growth (AWM).
- Predicts community OD using:
  - Bayesian AWM × g (posterior model integration)
  - observed AWM × g (oracle ceiling)
  - equal-abundance baseline
  - trait-shuffled baseline
- Computes composition and OD performance metrics (RMSE, R², JS divergence, LOO where available).
- Saves posterior prediction objects for reuse.

#### Part A figure outputs (`results/figures/`)
- `Fig_3a.tiff` — single-panel composition fit (predicted vs observed relative abundance)
- `Fig_3b.tiff` — OD predictions using Bayesian AWM × g
- `Fig_S4.tiff` — composition fits facetted by stress
- `Fig_S5a.tiff` — OD predictions using observed AWM × g (oracle)
- `Fig_S5b.tiff` — OD predictions using equal-abundance baseline
- `Fig_S5c.tiff` — OD predictions using trait-shuffled baseline

#### Part A table outputs (`results/tables/`)
- `Table_S9_composition_fit_metrics_by_stress.csv`
- `Table_S10_loo_elpd_by_stress.csv` *(if `log_lik` / LOO is available)*
- `Table_S11_model_overall_metrics.csv`
- `Table_S12_js_by_stress.csv`
- `Table_S13_loo_summary_stats.csv` *(if LOO is available)*
- `Table_S14_stan_keyparam_diagnostics.csv`
- `Table_SXX_stan_sampler_diagnostics.csv` *(placeholder numbering in current script)*
- `Table_SXX_composition_interval_coverage90_by_stress.csv` *(placeholder numbering in current script)*

#### Part A RDS outputs (`results/rds/`)
- `bayes_pred_comp.rds`
- `bayes_awms.rds`

#### Part B (blocked cross-validation) overview
- Evaluates composition and OD prediction performance on held-out test folds.
- Supports two blocking strategies:
  - `BLOCK_TYPE = "Com_Id"` (blocked K-fold by community ID; default/recommended)
  - `BLOCK_TYPE = "Stress"` (leave-one-stress-out)
- Writes summary tables with a dynamic filename tag based on block type and whether taxon biases are enabled.

#### Part B figure outputs (`results/figures/`)
- `Fig_S6.tiff` — CV test composition fits by stress
- `Fig_S7.tiff` — CV test OD predictions

#### Part B table outputs (`results/tables/`)
Filenames include a generated tag:
`blockedcv_<block_type>_bias<0/1>`

Examples (actual names depend on settings):
- `Table_S15_<tag>_by_fold.csv`
- `Table_S16_<tag>_by_stress.csv`
- `Table_S17_<tag>_od_by_fold_summary.csv`
- `Table_S18_<tag>_summary.csv`

---

## Running the pipeline

There is **no `analysis.R` or `run_all.R` in the current script set you shared**. Run scripts manually in order from the project root.

Recommended order:

```r
# 1) Preprocess ASV / metadata objects
source("scripts/01_data_preprocessing.R")

# 2) Fig 1 panels + tables
source("scripts/02_fig1_diversity_permanova.R")
source("scripts/03_fig1_indicator_species.R")

# 3) Fig 2 growth z-scores
source("scripts/04_fig2_growth_zscores.R")

# 4) Fig 3 Bayesian composition + OD + blocked CV
source("scripts/05_fig3_bayes_od.R")
```

If your scripts are in the project root (not `scripts/`), either move them into `scripts/` or adjust the `source()` calls accordingly.

