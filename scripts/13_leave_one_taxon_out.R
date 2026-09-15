# ======================================================================
# 13_leave_one_taxon_out.R
#
# KILL TEST 2: does the growth-to-abundance mapping predict taxa whose
# community behaviour was never used for calibration?
#
# For each taxon j: the model is fitted only on communities that do NOT
# contain j, then predicts the communities that do contain j. Taxon j's
# bias term delta_j is therefore informed only by its shrinkage prior
# (it collapses towards 0), so its predicted abundance rests on its
# monoculture growth rate alone. (For the reference taxon, whose delta
# is fixed at 0 by the identifiability constraint, the held-out
# prediction is structurally identical.)
#
# Success (held-out taxon abundances track observations) upgrades the
# claim from "interpolation among familiar taxa" to "a transferable
# trait rule". Failure establishes that taxon-specific calibration is
# required.
#
# Outputs:
#   results/tables/Table_S27_leave_one_taxon_out.csv
#   results/rds/bayes_loto_<taxon>.rds  (12 fits, cached)
# ======================================================================

source("scripts/utils_bayes_prep.R")

CHAINS <- 4; ITER <- 2000; WARMUP <- 1000
ADAPT  <- list(adapt_delta = 0.99, max_treedepth = 15)
CORES  <- 2
NDRAWS <- 250

message("[1/3] Preparing data ...")
trait_r <- growth_table_default()
prep    <- prep_bayes_data(trait_r)
design  <- design_communities()

message("[2/3] Compiling Stan model ...")
sm <- compile_softmax("softmax_dirichlet_refit.stan")

message("[3/3] Leave-one-taxon-out over ", length(prep$taxa_levels), " taxa ...")

rows      <- vector("list", length(prep$taxa_levels))
focal_all <- vector("list", length(prep$taxa_levels))

for (j in seq_along(prep$taxa_levels)) {
  tx <- prep$taxa_levels[j]

  coms_with    <- unique(design$Com_Id[design$Id == tx])
  train_ids <- prep$SAMPLES$SampleID[!(prep$SAMPLES$Com_Id %in% coms_with)]
  test_ids  <- prep$SAMPLES$SampleID[  prep$SAMPLES$Com_Id %in% coms_with]

  message("Taxon ", tx, " | train samples = ", length(train_ids),
          " | test samples = ", length(test_ids))
  if (length(train_ids) < 20) { message("  too few training samples, skipping"); next }

  train_obj <- build_ragged(prep$ord_full, train_ids, prep$taxa_levels,
                            prep$stress_lvls, prep$rich_lvls, FALSE, 1L)
  test_obj  <- build_ragged(prep$ord_full, test_ids,  prep$taxa_levels,
                            prep$stress_lvls, prep$rich_lvls, FALSE, 1L)

  fit <- fit_softmax_cached(
    sm, train_obj$stan_data,
    cache_path = P_RDS(sprintf("bayes_loto_%s.rds", tx)),
    chains = CHAINS, iter = ITER, warmup = WARMUP,
    adapt = ADAPT, seed = 4000 + j, cores = CORES
  )
  draws <- rstan::extract(fit)
  nd <- dim(draws$kappa)[1]
  set.seed(5000 + j)
  sel <- if (nd > NDRAWS) sample.int(nd, NDRAWS) else seq_len(nd)

  pred <- predict_comp_and_awm(test_obj, draws, sel, FALSE, prep$taxa_levels)$pred_comp

  focal <- pred %>% dplyr::filter(Id == tx)
  other <- pred %>% dplyr::filter(Id != tx)

  rows[[j]] <- tibble::tibble(
    taxon            = tx,
    n_train_com      = length(unique(design$Com_Id)) - length(coms_with),
    n_test_samples   = length(test_ids),
    focal_RMSE       = round(sqrt(mean((focal$p_obs - focal$p_hat)^2)), 3),
    focal_bias       = round(mean(focal$p_hat - focal$p_obs), 3),
    focal_spearman   = round(suppressWarnings(
                         cor(focal$p_obs, focal$p_hat, method = "spearman")), 3),
    others_RMSE      = round(sqrt(mean((other$p_obs - other$p_hat)^2)), 3)
  )
  focal_all[[j]] <- focal %>% dplyr::mutate(held_out_taxon = tx)
}

tab    <- dplyr::bind_rows(rows)
pooled <- dplyr::bind_rows(focal_all)

wm_focal <- weighted_comp_metrics(pooled)
overall <- tibble::tibble(
  taxon          = "POOLED (all held-out taxa)",
  n_train_com    = NA, n_test_samples = nrow(pooled),
  focal_RMSE     = round(sqrt(mean((pooled$p_obs - pooled$p_hat)^2)), 3),
  focal_bias     = round(mean(pooled$p_hat - pooled$p_obs), 3),
  focal_spearman = round(suppressWarnings(
                     cor(pooled$p_obs, pooled$p_hat, method = "spearman")), 3),
  others_RMSE    = NA
)
tab <- dplyr::bind_rows(tab, overall)

readr::write_csv(tab, P_TAB("Table_S27_leave_one_taxon_out.csv"))
saveRDS(pooled, P_RDS("loto_focal_predictions.rds"))
print(as.data.frame(tab))
message("Pooled focal weighted R2 = ", round(wm_focal["wR2"], 3),
        " | Pearson r = ", round(cor(pooled$p_obs, pooled$p_hat), 3))
message("Done: Table S27 written.")
