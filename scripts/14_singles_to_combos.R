# ======================================================================
# 14_singles_to_combos.R
#
# KILL TEST 3: can community outcomes under COMBINED stressors be
# predicted from communities observed only under the control and
# SINGLE stressors?
#
# The composition model is trained on the Control, pH, Sal and Temp
# regimes only; the four combination regimes (pHSal, pHTemp, SalTemp,
# pHSalTemp) are predicted with their stress-specific slopes shrunk to
# the hierarchical mean (no community data from any combination regime
# enters the fit). Monoculture growth rates measured under the
# combination regimes are used as trait inputs, as in the LOSO analysis.
#
# The OD regression is trained on the single-stressor regimes without
# stress terms (OD ~ AWM + richness) and predicts combination-regime OD.
#
# Success justifies a genuine multiple-stressor forecasting claim
# (combination outcomes from single-stressor community calibration).
#
# Outputs:
#   results/tables/Table_S28_singles_to_combos.csv
#   results/rds/bayes_singles_fit.rds  (1 fit, cached)
# ======================================================================

source("scripts/utils_bayes_prep.R")

CHAINS <- 4; ITER <- 2000; WARMUP <- 1000
ADAPT  <- list(adapt_delta = 0.99, max_treedepth = 15)
CORES  <- 2
NDRAWS <- 250

SINGLES <- c("Control", "pH", "Sal", "Temp")

message("[1/3] Preparing data ...")
trait_r <- growth_table_default()
prep    <- prep_bayes_data(trait_r)
comm_od <- load_comm_od()

combos    <- setdiff(prep$stress_lvls, SINGLES)
train_ids <- prep$SAMPLES$SampleID[prep$SAMPLES$Stress %in% SINGLES]
test_ids  <- prep$SAMPLES$SampleID[prep$SAMPLES$Stress %in% combos]
message("Train (singles): ", length(train_ids), " samples | Test (combos): ",
        length(test_ids), " samples")

message("[2/3] Compiling and fitting on single-stressor regimes ...")
sm <- compile_softmax("softmax_dirichlet_refit.stan")

train_obj <- build_ragged(prep$ord_full, train_ids, prep$taxa_levels,
                          prep$stress_lvls, prep$rich_lvls, FALSE, 1L)
test_obj  <- build_ragged(prep$ord_full, test_ids,  prep$taxa_levels,
                          prep$stress_lvls, prep$rich_lvls, FALSE, 1L)

fit <- fit_softmax_cached(
  sm, train_obj$stan_data,
  cache_path = P_RDS("bayes_singles_fit.rds"),
  chains = CHAINS, iter = ITER, warmup = WARMUP,
  adapt = ADAPT, seed = 6000, cores = CORES
)
draws <- rstan::extract(fit)
nd <- dim(draws$kappa)[1]
set.seed(6001)
sel <- if (nd > NDRAWS) sample.int(nd, NDRAWS) else seq_len(nd)

message("[3/3] Predicting combination regimes ...")
pred_test  <- predict_comp_and_awm(test_obj,  draws, sel, FALSE, prep$taxa_levels)
pred_train <- predict_comp_and_awm(train_obj, draws, sel, FALSE, prep$taxa_levels)

# --- Composition accuracy per combination regime ---
per_sample <- pred_test$pred_comp %>%
  dplyr::group_by(SampleID, Stress) %>%
  dplyr::summarise(
    RMSE = sqrt(mean((p_obs - p_hat)^2)),
    JS   = jsd_local(p_obs, p_hat),
    dom  = Id[which.max(p_hat)] == Id[which.max(p_obs)],
    .groups = "drop"
  )

by_regime <- per_sample %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(n = dplyr::n(), comp_RMSE = mean(RMSE), comp_JS = mean(JS),
                   dom_accuracy = mean(dom), .groups = "drop")

wm <- weighted_comp_metrics(pred_test$pred_comp)

# --- OD: train on singles (no stress terms), predict combos ---
agg_awm <- function(x) x %>%
  dplyr::group_by(Com_Id, Stress, Diversity) %>%
  dplyr::summarise(AWM = mean(AWM_mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Stress = as.character(Stress), Diversity = as.character(Diversity))

od_train <- comm_od %>%
  dplyr::inner_join(agg_awm(pred_train$awm_post), by = c("Com_Id","Stress","Diversity")) %>%
  dplyr::mutate(Diversity = factor(Diversity, levels = as.character(prep$rich_lvls)))
od_test <- comm_od %>%
  dplyr::inner_join(agg_awm(pred_test$awm_post), by = c("Com_Id","Stress","Diversity")) %>%
  dplyr::mutate(Diversity = factor(Diversity, levels = as.character(prep$rich_lvls)))

m_od    <- lm(OD ~ AWM + Diversity, data = od_train)
pred_od <- as.numeric(predict(m_od, newdata = od_test))
od_r2   <- 1 - sum((od_test$OD - pred_od)^2) / sum((od_test$OD - mean(od_test$OD))^2)
od_rmse <- sqrt(mean((od_test$OD - pred_od)^2))

summary_row <- tibble::tibble(
  Stress = "POOLED (all combos)", n = nrow(per_sample),
  comp_RMSE = mean(per_sample$RMSE), comp_JS = mean(per_sample$JS),
  dom_accuracy = mean(per_sample$dom)
)

tab <- dplyr::bind_rows(by_regime, summary_row) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::mutate(comp_wR2_pooled = c(rep(NA, nrow(by_regime)), round(wm["wR2"], 3)),
                OD_R2_pooled    = c(rep(NA, nrow(by_regime)), round(od_r2, 3)),
                OD_RMSE_pooled  = c(rep(NA, nrow(by_regime)), round(od_rmse, 3)))

readr::write_csv(tab, P_TAB("Table_S28_singles_to_combos.csv"))
print(as.data.frame(tab))
message("Pooled combo composition: wR2 = ", round(wm["wR2"], 3),
        ", wRMSE = ", round(wm["wRMSE"], 3))
message("Done: Table S28 written.")
