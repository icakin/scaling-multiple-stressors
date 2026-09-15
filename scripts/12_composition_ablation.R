# ======================================================================
# 12_composition_ablation.R
#
# KILL TEST 1: does growth add out-of-sample information beyond taxon
# identity? Ablates the composition model under blocked community-
# identity cross-validation:
#
#   full        : eta = kappa[s] * g_z + delta[taxon]   (as published)
#   growth_only : eta = kappa[s] * g_z                  (taxon biases off)
#   biases_only : eta = delta[taxon]                    (growth zeroed)
#
# All three variants are fitted with the same Stan model: growth_only
# sets use_taxon_biases = 0; biases_only sets g_z = 0 in the data (the
# kappa parameters are then unidentified and sample from their prior,
# which is harmless but may produce low-ESS warnings for kappa; ignore
# those for that variant).
#
# Evaluation on held-out communities (same 5-fold Com_Id blocking as
# script 05 Part B): pooled out-of-fold weighted R2/RMSE, mean per-sample
# RMSE and Jensen-Shannon divergence, dominant-taxon accuracy, and the
# pooled held-out log predictive density (lpd; higher is better), which
# is the cross-validated analogue of ELPD.
#
# Interpretation: if biases_only approaches full, taxon identity (not
# growth) carries the composition prediction. If growth_only is strong
# and full clearly beats biases_only, growth carries transferable
# information.
#
# Outputs:
#   results/tables/Table_S26_composition_ablation.csv
#   results/rds/bayes_ablation_<variant>_fold<k>.rds  (15 fits, cached)
# ======================================================================

source("scripts/utils_bayes_prep.R")

CHAINS <- 4; ITER <- 2000; WARMUP <- 1000
ADAPT  <- list(adapt_delta = 0.99, max_treedepth = 15)
CORES  <- 2
NDRAWS <- 250
K_FOLDS <- 5

# ---------------------- Data & folds ---------------------------------
message("[1/3] Preparing data ...")
trait_r <- growth_table_default()
prep    <- prep_bayes_data(trait_r)

# Same blocked fold construction as script 05 Part B
set.seed(123)
blocks   <- sample(unique(prep$SAMPLES$Com_Id))
fold_map <- setNames(rep_len(seq_len(K_FOLDS), length(blocks)), blocks)
prep$SAMPLES$fold <- as.integer(fold_map[prep$SAMPLES$Com_Id])

message("[2/3] Compiling Stan model ...")
sm <- compile_softmax("softmax_dirichlet_refit.stan")

# ---------------------- Helpers --------------------------------------
dirichlet_lpdf <- function(x, alpha) {
  lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
}
log_mean_exp <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }

# Held-out predictions + log predictive density for one fitted fold
evaluate_test <- function(obj, draws, sel, taxa_levels) {
  ord   <- obj$ord_sub
  smp   <- obj$samples_sub
  g_z_v <- as.numeric(ord$g_z)
  out_pred <- vector("list", nrow(smp))
  lpd      <- numeric(nrow(smp))
  for (n in seq_len(nrow(smp))) {
    seg     <- obj$idx_split[[n]]
    s       <- obj$stress_id[n]
    taxa_ix <- obj$tax_of[seg]
    pseg    <- ord$Abundance[seg]
    pj      <- pseg + 1e-12; pj <- pj / sum(pj)   # same jitter as the model
    P   <- matrix(NA_real_, length(sel), length(seg))
    lp  <- numeric(length(sel))
    for (ii in seq_along(sel)) {
      d      <- sel[ii]
      kap    <- get_kappa_cv(draws, d, s, 1L)
      p      <- softmax_vec(kap * g_z_v[seg] + draws$delta[d, taxa_ix])
      P[ii,] <- p
      lp[ii] <- dirichlet_lpdf(pj, draws$phi[d, s] * p)
    }
    lpd[n] <- log_mean_exp(lp)
    out_pred[[n]] <- tibble::tibble(
      SampleID = smp$SampleID[n],
      Id       = taxa_levels[taxa_ix],
      p_obs    = pseg,
      p_hat    = colMeans(P)
    )
  }
  list(pred = dplyr::bind_rows(out_pred), lpd = lpd)
}

# ---------------------- Ablation loop --------------------------------
message("[3/3] Running ablation x blocked CV ...")

variants <- list(
  full        = list(zero_g = FALSE, biases = 1L),
  growth_only = list(zero_g = FALSE, biases = 0L),
  biases_only = list(zero_g = TRUE,  biases = 1L)
)

rows <- vector("list", length(variants))

for (vn in names(variants)) {
  v <- variants[[vn]]
  message("Variant: ", vn)

  ord_v <- prep$ord_full
  if (v$zero_g) ord_v$g_z <- 0

  pred_all <- vector("list", K_FOLDS)
  lpd_all  <- c()

  for (k in seq_len(K_FOLDS)) {
    test_ids  <- prep$SAMPLES$SampleID[prep$SAMPLES$fold == k]
    train_ids <- prep$SAMPLES$SampleID[prep$SAMPLES$fold != k]

    train_obj <- build_ragged(ord_v, train_ids, prep$taxa_levels,
                              prep$stress_lvls, prep$rich_lvls, FALSE, v$biases)
    test_obj  <- build_ragged(ord_v, test_ids,  prep$taxa_levels,
                              prep$stress_lvls, prep$rich_lvls, FALSE, v$biases)

    fit <- fit_softmax_cached(
      sm, train_obj$stan_data,
      cache_path = P_RDS(sprintf("bayes_ablation_%s_fold%d.rds", vn, k)),
      chains = CHAINS, iter = ITER, warmup = WARMUP,
      adapt = ADAPT, seed = 900 + k, cores = CORES
    )
    draws <- rstan::extract(fit)
    nd <- dim(draws$kappa)[1]
    set.seed(3000 + k)
    sel <- if (nd > NDRAWS) sample.int(nd, NDRAWS) else seq_len(nd)

    ev <- evaluate_test(test_obj, draws, sel, prep$taxa_levels)
    pred_all[[k]] <- ev$pred
    lpd_all       <- c(lpd_all, ev$lpd)
  }

  pooled <- dplyr::bind_rows(pred_all)
  wm     <- weighted_comp_metrics(pooled)

  per_sample <- pooled %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      RMSE = sqrt(mean((p_obs - p_hat)^2)),
      JS   = jsd_local(p_obs, p_hat),
      dom  = Id[which.max(p_hat)] == Id[which.max(p_obs)],
      .groups = "drop"
    )

  rows[[vn]] <- tibble::tibble(
    variant       = vn,
    comp_wR2_oof  = round(wm["wR2"], 3),
    comp_wRMSE_oof= round(wm["wRMSE"], 3),
    RMSE_mean     = round(mean(per_sample$RMSE), 3),
    JS_mean       = round(mean(per_sample$JS), 3),
    dom_accuracy  = round(mean(per_sample$dom), 3),
    lpd_total     = round(sum(lpd_all), 1),
    lpd_mean      = round(mean(lpd_all), 3),
    n_test        = nrow(per_sample)
  )
}

tab <- dplyr::bind_rows(rows)
readr::write_csv(tab, P_TAB("Table_S26_composition_ablation.csv"))
print(as.data.frame(tab))
message("Done: Table S26 written.")
