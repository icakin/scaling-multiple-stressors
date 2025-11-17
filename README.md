# Scaling the impacts of multiple stressors from populations to ecosystems

This repository contains the R scripts used in the manuscript:  
Carmichael H., Cakin I., Busi S.B., Read D., Yvon-Durocher G. “Scaling the impacts of multiple stressors from populations to ecosystems”.

## Contents
- `analysis.R` — runs all analyses described in the manuscript, including diversity–biomass models, db-RDA, PERMANOVA, indicator species analysis, Shannon diversity, trait processing, and the Bayesian softmax (trait → composition) model with biomass prediction.  
- (optional) `figures/`, `tables/`, `data/` — folders for outputs and input data.

## R version and packages
- R ≥ 4.3  
- Required packages: tidyverse, vegan, phyloseq, indicspecies, ggplot2, ggpubr, gridExtra, reshape2, broom, rstatix, car, emmeans, AICcmodavg, lmtest, magrittr, scales, ggrepel, patchwork, ape, pairwiseAdonis, ggtree, stringr, picante, brms, rstan.

## Reproducible setup
Install packages, then run:

```r
# from R
# install.packages("renv")
renv::init()
renv::restore()
source("analysis.R")
```

---

## Analysis pipeline

### 1. Data preprocessing

**Script:** `scripts/01_data_preprocessing.R`  
**Outputs:**

- `results/rds/phyloseq_filtered.rds`
- `results/rds/ASV_matrix_clean.rds`
- `results/rds/metadata_clean.rds`
- `results/rds/taxa_table.rds`

**Summary:**  
Reads in ASV counts, taxonomy, sample metadata, and the phylogenetic tree, builds a `phyloseq` object, converts to relative abundance, filters low-abundance ASVs, subsets to Eight-species, non-evolved communities, and exports aligned abundance and metadata matrices for all downstream scripts.

---

### 2. Fig 1 – Diversity, PERMANOVA & db-RDA

**Script:** `scripts/02_fig1_diversity_permanova.R`  
**Main products:**

- **Figures**
  - `results/figures/Fig_1a_Shannon_diversity_boxplot_with_signif.tiff`
  - `results/figures/Fig_1b_dbRDA_gradients.tiff`
  - `results/figures/Fig_1c_PERMANOVA_pairwise_heatmap.tiff`
- **Tables**
  - `results/tables/Table_S4_Tukey_posthoc_effectsizes.csv`
  - `results/tables/Table_S5_PERMANOVA_pairwise.csv`
  - `results/tables/Shannon_ANOVA_effectsizes.csv`

**Summary:**  

- Computes Bray–Curtis distances from the cleaned ASV matrix.  
- Runs db-RDA (distance-based RDA) against pH, salinity, and temperature and plots stress-wise centroids/ellipses (Fig 1b).  
- Performs global PERMANOVA and pairwise PERMANOVA across stress treatments; visualises R² values in a heatmap (Fig 1c, Table S5).  
- Calculates Shannon diversity, fits ANOVA + Tukey post-hoc tests, and exports effect sizes (Fig 1a, Table S4, Shannon_ANOVA_effectsizes.csv).

---

### 3. Fig 1 – Indicator species

**Script:** `scripts/03_fig1_indicator_species.R`  
**Main products:**

- **Figure**
  - `results/figures/Fig_1d_ISA_Top1_1vsRest_singlepanel.tiff`
- **Table**
  - `results/tables/Table_S6_ISA_stats_full.csv`

**Summary:**  

- Uses indicator species analysis (`indicspecies::multipatt`, IndVal.g) to identify taxa most associated with each stress treatment.  
- Cleans and harmonises genus-level taxonomy (including disambiguating Pseudomonas ASVs).  
- Extracts the top indicator taxon per stress and plots the indicator value by genus and condition (Fig 1d).  
- Writes the full indicator statistics for all ASVs × stress groups (Table S6).

---

### 4. Fig 2 – Growth z-scores across stressors

**Script:** `scripts/04_fig2_growth_zscores.R`  
**Main products:**

- **Figure**
  - `results/figures/Fig_2_DOTPLOT_z_by_stress_compact.tiff`
- **Tables**
  - `results/tables/Table_S7_Growth_all_sources_stressors_long.csv`
  - `results/tables/Table_S8_Zscores_summary.csv`

**Summary:**  

- Loads single-taxon growth-rate fits from multiple experiments (Temp, pH, Salinity).  
- Harmonises columns and normalises experimental levels to discrete stress categories (Control, Temp, pH, Sal, and their combinations).  
- Combines all sources into a long-format growth dataset and filters non-negative rates (Table S7).  
- Computes within-stress z-scores for each taxon and summarises mean ± SE across replicates (Table S8).  
- Produces a compact multi-panel dotplot of z-scored growth rates by taxon and stress (Fig 2).

---

### 5. Fig 3 – Bayesian composition model & OD predictions

**Script:** `scripts/05_fig3_bayes_od.R`  
**Main products:**

- **Figures**
  - `results/figures/Fig_3a.tiff`  
  - `results/figures/Fig_3b.tiff`  
  - `results/figures/Fig_S1.tiff`  
  - `results/figures/Fig_S2a.tiff`  
  - `results/figures/Fig_S2b.tiff`  
  - `results/figures/Fig_S2c.tiff`
- **Tables**
  - `results/tables/Table_S9_OD_metrics_by_stress.csv`
  - `results/tables/Table_S10_loo_elpd_by_stress.csv`
  - `results/tables/Table_S11_model_overall_metrics.csv`
  - `results/tables/Table_S12_js_by_stress.csv`
  - `results/tables/Table_S13_loo_summary_stats.csv`
- **RDS**
  - `results/rds/bayes_pred_comp.rds`
  - `results/rds/bayes_awms.rds`

**Summary:**  

- Constructs a Dirichlet–softmax Bayesian model linking taxon-level z-scored growth to observed relative abundances.  
- Fits with `rstan`, generates posterior predicted compositions, and calculates abundance-weighted mean growth (AWM).  
- Loads community OD data and tests four predictors of OD (observed AWM×g, Bayesian AWM×g, equal-abundance baseline, random baseline).  
- Summarises model performance using RMSE, JS distance, and LOO metrics.

---

## Running the pipeline

From the project root in R:

```r
source("run_all.R")
```
