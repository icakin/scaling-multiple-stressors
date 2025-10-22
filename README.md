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
# renv::init(); renv::restore()
source("analysis.R")
