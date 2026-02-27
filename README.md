# Stress-specific growth rates predict microbial community assembly and biomass yield under concurrent environmental stressors

This repository contains all code, data, and manuscript source files for:

**Carmichael H.\*, Cakin I.\*, Busi S.B., Read D. & Yvon-Durocher G.**
*Stress-specific growth rates predict microbial community assembly and biomass yield under concurrent environmental stressors.*

\* These authors contributed equally.

## Overview

We show that stress-specific monoculture growth rates predict how synthetic bacterial communities assemble and how much biomass they produce under factorial combinations of temperature, pH and salinity. A Bayesian softmax model maps taxon-level growth rates to relative abundances across stress regimes, and an abundance-weighted mean growth metric links predicted composition to endpoint biomass yield (OD at 600 nm). The analysis pipeline runs from raw amplicon-sequencing and growth-curve data through to publication-ready figures and a typeset manuscript.

## Repository structure

```
scaling-multiple-stressors/
|
|-- data/                          # Raw input data (CSV, TSV, Newick tree)
|-- scripts/                       # R analysis scripts (numbered pipeline)
|-- results/
|   |-- figures/                   # Output figures (PNG, TIFF)
|   |-- tables/                    # Output tables (CSV)
|   |-- rds/                       # Cached R objects (RDS)
|
|-- manuscript/                    # Quarto manuscript project
|   |-- manuscript.qmd             # Main text
|   |-- supplementary.qmd          # Supplementary material
|   |-- references.bib             # Bibliography
|   |-- nature-microbiology.csl    # Citation style
|   |-- header.tex                 # LaTeX preamble (main)
|   |-- header-supp.tex            # LaTeX preamble (supplementary)
|   |-- _quarto.yml                # Quarto project configuration
|   |-- _output/                   # Rendered PDFs and DOCX (not tracked)
|
|-- softmax_dirichlet.stan         # Bayesian composition model (Stan)
|-- softmax_dirichlet_blockedcv.stan  # Same model, cross-validation variant
|-- run_all.R                      # Entry point: runs the full pipeline
|-- Rproject.Rproj                 # RStudio project file
|-- LICENSE                        # MIT License
```

## Prerequisites

### Software

- **R** >= 4.3 (tested with R 4.4)
- **RStan** and a working C++ toolchain (for Bayesian model fitting)
- **Quarto** >= 1.4 (for manuscript rendering)
- **TeX Live** or equivalent LaTeX distribution with LuaLaTeX (for PDF output)

### R packages

All required packages are installed automatically the first time you run any script, via the `ensure_packages()` function in `scripts/utils_functions.R`. The main dependencies are:

| Category | Packages |
|----------|----------|
| Data wrangling | `dplyr`, `tidyr`, `tibble`, `readr`, `stringr`, `reshape2`, `purrr` |
| Plotting | `ggplot2`, `ggpubr`, `ggrepel`, `patchwork`, `scales`, `ragg`, `ggtree` |
| Microbial ecology | `phyloseq`, `vegan`, `pairwiseAdonis`, `indicspecies`, `microeco`, `ape`, `picante`, `phangorn` |
| Statistics | `broom`, `car`, `effectsize`, `emmeans`, `permute` |
| Bayesian modelling | `rstan`, `loo`, `philentropy` |
| Growth curves | `nls.multstart`, `MuMIn` |

## Input data

All raw data files live in `data/`:

| File | Description |
|------|-------------|
| `ASV.counts.tsv` | ASV count matrix (samples x taxa) |
| `ASV.tax.tsv` | Taxonomic assignments for each ASV |
| `Sequencing_metadata2.csv` | Sample metadata (stress treatment, richness, community identity) |
| `ASV1-13_tree.txt` | Phylogenetic tree (Newick format) |
| `ALl_Abundance_Data.csv` | Observed community relative abundances |
| `CommunityOD_All_Jan23_edit_zeros.csv` | Community-level OD measurements |
| `LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv` | Monoculture growth-curve fits (pH gradient) |
| `LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv` | Monoculture growth-curve fits (temperature gradient) |
| `LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv` | Monoculture growth-curve fits (salinity gradient) |
| `Cut_OD_data/` | Per-taxon OD time-series for growth-curve fitting |

See `data/README_DATA.txt` for further details on data provenance and column definitions.

## Analysis pipeline

### Quick start

From the project root in R (or RStudio with `Rproject.Rproj` open):

```r
source("run_all.R")
```

This executes scripts 01-05 in order and saves a `results/sessionInfo.txt` log. The full pipeline takes approximately 30-60 minutes depending on hardware, with the Bayesian model fitting in script 05 accounting for most of the runtime. Stan model fits are cached as RDS files; subsequent runs skip refitting if the cache exists.

### Step-by-step

Each script is self-contained: it loads its own dependencies via `source("scripts/utils_functions.R")` and reads inputs from `data/` or `results/rds/`. Scripts must be run in numerical order because later scripts depend on RDS objects produced by earlier ones.

#### Script 01 -- Data preprocessing

**File:** `scripts/01_data_preprocessing.R`

Loads raw ASV counts, taxonomy, sample metadata, and the phylogenetic tree; constructs a `phyloseq` object; converts to relative abundance; filters taxa below 0.5% mean abundance; optionally merges two *Pseudomonas* ASVs; and subsets to eight-species, non-evolved communities. Outputs four aligned RDS objects used by all downstream scripts.

**Outputs:** `results/rds/{phyloseq_filtered, ASV_matrix_clean, metadata_clean, taxa_table}.rds`

#### Script 02 -- Diversity, db-RDA and PERMANOVA (Fig. 1a-c)

**File:** `scripts/02_fig1_diversity_permanova.R`

Computes Shannon diversity across stress treatments (ANOVA + Tukey post-hoc with effect sizes), fits a distance-based redundancy analysis (db-RDA) on Bray-Curtis distances using pH, salinity and temperature as constraints, and runs global and pairwise PERMANOVA. Produces three figure panels and associated supplementary tables. Also caches plot objects as RDS for composite figure assembly in script 03.

**Outputs:**
- Figures: `Fig_1a_Shannon_diversity_boxplot_with_signif`, `Fig_1b_dbRDA_gradients`, `Fig_1c_PERMANOVA_pairwise_heatmap`
- Tables: `Table_S5_Tukey_posthoc_effectsizes.csv`, `Table_S6_PERMANOVA_pairwise.csv`, `Shannon_ANOVA_effectsizes.csv`
- Cached plots: `results/rds/{fig1a_plot, fig1b_plot, fig1c_plot}.rds`

#### Script 03 -- Indicator species analysis (Fig. 1d) and composite Fig. 1

**File:** `scripts/03_fig1_indicator_species.R`

Runs indicator value analysis (IndVal.g via `indicspecies::multipatt`) with one-vs-rest contrasts for each of the eight stress conditions. Extracts the top indicator genus per condition. Loads cached panels a-c from script 02 and assembles the four-panel composite Fig. 1 using `patchwork`.

**Outputs:**
- Figures: `Fig_1d_ISA_Top1_1vsRest_singlepanel`, `Fig_1` (composite)
- Table: `Table_S7_ISA_stats_full.csv`

#### Script 04 -- Growth-rate z-scores (Fig. 2)

**File:** `scripts/04_fig2_growth_zscores.R`

Loads monoculture growth-rate fits from three stressor gradients (temperature, pH, salinity), harmonises column names and condition labels, maps each combination to one of the eight factorial stress regimes, computes within-stress z-scores (number of standard deviations from the treatment mean), and produces the compact dotplot (Fig. 2). Also exports the unified long-format growth table and z-score summary used by script 05.

**Outputs:**
- Figure: `Fig_2_DOTPLOT_z_by_stress_compact`
- Tables: `Table_S8a_Growth_all_sources_stressors_long.csv`, `Table_S8_Zscores_summary.csv`
- RDS: `Growth_all_sources_stressors_long.rds`, `Zscores_objects.rds`

#### Script 05 -- Bayesian composition model, OD prediction and cross-validation (Fig. 3)

**File:** `scripts/05_fig3_bayes_od.R`

This is the most complex script, split into two parts:

**Part A -- Model fitting and prediction.** Fits a Dirichlet-softmax Bayesian model (`softmax_dirichlet.stan`) via RStan that maps stress-standardised growth rates to community composition. The model estimates stress-regime-specific scaling parameters (kappa), optional taxon biases (delta), and a concentration parameter (phi). Posterior predicted compositions are converted to abundance-weighted mean growth (AWM), which is then regressed against observed OD to predict community biomass. Predictions are compared against three baselines: equal-abundance, random-Dirichlet, and an oracle using observed compositions. Stan diagnostics (divergences, R-hat, ESS) are exported.

**Part B -- Blocked cross-validation.** Performs five-fold CV blocked by community identity (Com_Id): in each fold the composition model and OD regression are fit on training communities only, and composition accuracy (RMSE, Jensen-Shannon divergence) and OD accuracy (RMSE, R-squared) are evaluated on held-out communities. Stress-specific standardisation parameters are computed from training data only to prevent information leakage.

**Outputs:**
- Figures: `Fig_3a`, `Fig_3b`, `Fig_3` (composite), `Fig_S4` through `Fig_S7`
- Tables: `Table_S9` through `Table_S18` (composition metrics, LOO, JS divergence, model comparison, Stan diagnostics, CV summaries)
- RDS: `bayes_pred_comp.rds`, `bayes_awms.rds`, `bayes_fit_partA.rds`, `bayes_cv_fold[1-5].rds`

### Supplementary growth-curve scripts

These two scripts support the monoculture growth-curve fitting that produced the input data for the main pipeline. They are not called by `run_all.R` but can be run independently.

#### Script 06 -- Growth model library

**File:** `scripts/06_growth_curve_models.R`

Defines parametric growth-curve functions (Baranyi, Gompertz, Buchanan, logistic -- each with and without lag) and Sharpe-Schoolfield temperature-response helpers. These functions are sourced by script 07.

#### Script 07 -- Growth-curve fitting

**File:** `scripts/07_growth_curves.R`

Reads per-taxon OD time-series from `data/Cut_OD_data/`, fits all candidate growth models via `nls.multstart`, selects the best model per curve by AICc, and exports fitted parameters and diagnostic plots (Figs. S1-S3). Results feed into the growth-rate CSV files in `data/` used by the main pipeline.

## Stan models

Two Stan files in the project root implement the Dirichlet-softmax composition model:

- **`softmax_dirichlet.stan`** -- Full model with `generated quantities` block that computes per-observation log-likelihoods for LOO-CV via the `loo` package.
- **`softmax_dirichlet_blockedcv.stan`** -- Identical model structure but without `generated quantities`, used during blocked cross-validation (where LOO is not needed and omitting it speeds up sampling).

Both models parameterise relative abundance as:

```
alpha[j] = phi * softmax(kappa[stress] * g_z[j] + delta[j])
y ~ dirichlet(alpha)
```

where `g_z[j]` is the z-scored growth rate of taxon j under the focal stress regime, `kappa` scales the growth-to-abundance mapping per stress, `delta` captures taxon-specific biases, and `phi` controls the concentration (precision) of the Dirichlet likelihood.

## Manuscript

The manuscript is a [Quarto manuscript project](https://quarto.org/docs/manuscripts/) in `manuscript/`. It consists of two documents:

- **`manuscript.qmd`** -- Main text (title, abstract, introduction, methods, results, discussion, conclusions)
- **`supplementary.qmd`** -- Supplementary material (additional figures, tables S1-S18, supplementary methods)

Both documents share `references.bib` (BibTeX bibliography) and `nature-microbiology.csl` (Nature Microbiology citation style). Custom LaTeX preambles (`header.tex`, `header-supp.tex`) handle author affiliations, figure caption formatting ("**Fig. N |** Title"), line numbering, and supplementary figure numbering (S-prefix).

### Rendering the manuscript

From the `manuscript/` directory:

```bash
# Main manuscript (PDF and Word)
quarto render manuscript.qmd --to pdf
quarto render manuscript.qmd --to docx
```

The supplementary document must be rendered outside the Quarto manuscript project context (which bundles it as a notebook rather than a standalone PDF). A portable approach:

```bash
PROJDIR=$(pwd)/..
TMPDIR=$(mktemp -d)
mkdir -p "$TMPDIR/manuscript"
cp supplementary.qmd references.bib nature-microbiology.csl header-supp.tex "$TMPDIR/manuscript/"
ln -s "$PROJDIR/results" "$TMPDIR/results"
cd "$TMPDIR/manuscript"
quarto render supplementary.qmd --to pdf
quarto render supplementary.qmd --to docx
cp supplementary.pdf supplementary.docx "$PROJDIR/manuscript/_output/"
```

Rendered outputs are written to `manuscript/_output/` (not tracked by git).

## Shared utilities

`scripts/utils_functions.R` is sourced by every analysis script and provides:

- **`ensure_packages()`** -- Installs missing CRAN and Bioconductor packages on first run.
- **Path helpers** -- `P_IN()`, `P_TAB()`, `P_FIG()`, `P_RDS()` return paths into `data/`, `results/tables/`, `results/figures/`, and `results/rds/` respectively, creating directories if needed.
- **Stress labels and palettes** -- An eight-level factor (`Control`, `pH`, `Sal`, `Temp`, `pH x Sal`, `pH x Temp`, `Sal x Temp`, `pH x Sal x Temp`) with a consistent colour palette used across all figures.
- **Global settings** -- Random seed (123), Stan parallel cores, `stringsAsFactors = FALSE`.

## Reproducing the analysis

1. **Clone the repository** and open `Rproject.Rproj` in RStudio (or set the working directory to the project root).
2. **Run the pipeline:**
   ```r
   source("run_all.R")
   ```
   On first run this will install all required R packages, fit the Stan models (cached for subsequent runs), and produce all figures and tables in `results/`.
3. **Render the manuscript** (requires Quarto and a LaTeX distribution):
   ```bash
   cd manuscript
   quarto render manuscript.qmd --to pdf
   quarto render manuscript.qmd --to docx
   ```
4. **Render the supplementary** using the temp-directory approach described above.

## License

This project is released under the [MIT License](LICENSE).
