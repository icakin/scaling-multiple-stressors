# ======================================================================
# 09_loso_cv.R
#
# Leave-one-stress-out (LOSO) cross-validation of the trait-based
# composition -> biomass pipeline, as requested by Reviewer 4
# (Ecology Letters review, 2026): does the growth-rate mapping
# generalise across environments, not just across communities?
#
# Design: each of the 8 stress regimes is held out in turn. The
# Dirichlet-softmax model is fitted to the other 7 regimes only, and
# compositions are predicted for the held-out regime. Because the
# held-out regime contributes no data, its stress-specific slope
# kappa[s] and concentration phi[s] are informed only by the
# hierarchical priors (i.e. they shrink to the across-stress means),
# which is the honest prediction for an unseen environment.
#
# The OD regression is fitted on training regimes WITHOUT stress terms
# (OD ~ AWM + Diversity): stress-specific intercepts/slopes cannot be
# estimated for a regime absent from training, so the trait term has to
# carry the prediction. This makes LOSO a strictly harder test than the
# blocked Com_Id CV (Tables S15-S18).
#
# NOTE: within-stress z-scoring of growth rates uses only monoculture
# data (an exogenous input measured independently of the communities),
# so no community-level information leaks from the held-out regime.
#
# Outputs:
#   results/tables/Table_S20_losocv_by_stress.csv
#   results/tables/Table_S21_losocv_summary.csv
#   results/figures/Fig_S8_losocv.png / .tiff
#   results/rds/bayes_losocv_fold_<stress>.rds   (per-fold Stan fits, cached)
#
# Runtime: ~8 Stan fits (2000 iter, 4 chains each). Set CORES_CV to 4
# to run chains in parallel if your machine allows it.
# ======================================================================

source("scripts/utils_bayes_prep.R")
for (p in c("philentropy","ggplot2","ragg")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(ragg) })

# ---------------------- Options --------------------------------------
CHAINS_CV <- 4
ITER_CV   <- 2000
WARMUP_CV <- 1000
ADAPT_CV  <- list(adapt_delta = 0.99, max_treedepth = 15)
CORES_CV  <- 1          # set to 4 for parallel chains if stable on your machine
NDRAWS_COMP <- 250
NDRAWS_MI   <- 250

OD_FORMULA_LOSO <- stats::as.formula("OD ~ AWM + Diversity")

# ---------------------- Data -----------------------------------------
message("[1/3] Preparing data ...")
trait_r <- growth_table_default()
prep    <- prep_bayes_data(trait_r)
comm_od <- load_comm_od()

SAMPLES_cv  <- prep$SAMPLES
stress_lvls <- prep$stress_lvls
K_eff       <- length(stress_lvls)

message("[2/3] Compiling Stan model ...")
sm_cv <- compile_softmax(P_STAN("softmax_dirichlet_refit.stan"))

# ---------------------- LOSO loop ------------------------------------
message("[3/3] Running LOSO CV over ", K_eff, " stress regimes ...")

fold_rows   <- vector("list", K_eff)
od_pred_all <- vector("list", K_eff)

for (k in seq_len(K_eff)) {

  s_out    <- stress_lvls[k]
  test_ids  <- SAMPLES_cv$SampleID[SAMPLES_cv$Stress == s_out]
  train_ids <- SAMPLES_cv$SampleID[SAMPLES_cv$Stress != s_out]

  message("Fold ", k, "/", K_eff, " | held-out stress = ", s_out,
          " | train=", length(train_ids), " test=", length(test_ids))

  train_obj <- build_ragged(prep$ord_full, train_ids, prep$taxa_levels,
                            stress_lvls, prep$rich_lvls, FALSE, 1L)
  test_obj  <- build_ragged(prep$ord_full, test_ids,  prep$taxa_levels,
                            stress_lvls, prep$rich_lvls, FALSE, 1L)

  fit_cv <- fit_softmax_cached(
    sm_cv, train_obj$stan_data,
    cache_path = P_RDS(sprintf("bayes_losocv_fold_%s.rds", s_out)),
    chains = CHAINS_CV, iter = ITER_CV, warmup = WARMUP_CV,
    adapt = ADAPT_CV, seed = 700 + k, cores = CORES_CV
  )

  draws_cv <- rstan::extract(fit_cv)
  nd <- dim(draws_cv$kappa)[1]

  set.seed(1000 + k)
  sel_comp <- if (nd > NDRAWS_COMP) sample.int(nd, NDRAWS_COMP) else seq_len(nd)
  set.seed(2000 + k)
  sel_mi   <- if (nd > NDRAWS_MI)   sample.int(nd, NDRAWS_MI)   else seq_len(nd)

  # --- Composition accuracy on the held-out regime ---
  pred_test <- predict_comp_and_awm(test_obj, draws_cv, sel_comp, FALSE, prep$taxa_levels)

  comp_per_sample <- pred_test$pred_comp %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      RMSE = sqrt(mean((p_obs - p_hat)^2)),
      JS   = suppressWarnings(suppressMessages(philentropy::JSD(rbind(p_obs, p_hat)))),
      .groups = "drop"
    )

  # --- OD prediction: train regression on 7 regimes, predict held-out ---
  pred_train_mi <- predict_comp_and_awm(train_obj, draws_cv, sel_mi, FALSE, prep$taxa_levels)

  agg_awm <- function(x) x %>%
    dplyr::group_by(Com_Id, Stress, Diversity) %>%
    dplyr::summarise(AWM = mean(AWM_mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Stress = as.character(Stress), Diversity = as.character(Diversity))

  od_train <- comm_od %>%
    dplyr::inner_join(agg_awm(pred_train_mi$awm_post), by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::mutate(Diversity = factor(Diversity, levels = as.character(prep$rich_lvls)))

  od_test <- comm_od %>%
    dplyr::inner_join(agg_awm(pred_test$awm_post), by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::mutate(Diversity = factor(Diversity, levels = as.character(prep$rich_lvls)))

  od_rmse <- NA_real_; od_r2 <- NA_real_
  if (nrow(od_train) >= 10 && nrow(od_test) >= 3) {
    m_od    <- stats::lm(OD_FORMULA_LOSO, data = od_train)
    pred_od <- as.numeric(stats::predict(m_od, newdata = od_test))
    od_rmse <- sqrt(mean((od_test$OD - pred_od)^2))
    od_r2   <- 1 - sum((od_test$OD - pred_od)^2) / sum((od_test$OD - mean(od_test$OD))^2)
    od_pred_all[[k]] <- od_test %>% dplyr::mutate(pred_OD = pred_od, held_out = s_out)
  }

  fold_rows[[k]] <- tibble::tibble(
    held_out_stress = s_out,
    n_test          = length(test_ids),
    comp_RMSE       = mean(comp_per_sample$RMSE),
    comp_JS         = mean(comp_per_sample$JS),
    n_test_OD       = nrow(od_test),
    OD_RMSE         = od_rmse,
    OD_R2           = od_r2
  )
}

by_stress <- dplyr::bind_rows(fold_rows) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

summary_tbl <- by_stress %>%
  dplyr::summarise(
    comp_RMSE_mean = mean(comp_RMSE), comp_RMSE_sd = sd(comp_RMSE),
    comp_JS_mean   = mean(comp_JS),   comp_JS_sd   = sd(comp_JS),
    OD_RMSE_mean   = mean(OD_RMSE, na.rm = TRUE), OD_RMSE_sd = sd(OD_RMSE, na.rm = TRUE),
    OD_R2_mean     = mean(OD_R2,   na.rm = TRUE), OD_R2_sd   = sd(OD_R2,   na.rm = TRUE)
  ) %>%
  dplyr::mutate(Block_type = "Stress (LOSO)", K = K_eff, .before = 1) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

readr::write_csv(by_stress,   P_TAB("Table_S20_losocv_by_stress.csv"))
readr::write_csv(summary_tbl, P_TAB("Table_S21_losocv_summary.csv"))

# ---------------------- Figure S8 ------------------------------------
od_all <- dplyr::bind_rows(od_pred_all)

if (nrow(od_all) > 0) {
  rmse_pool <- sqrt(mean((od_all$OD - od_all$pred_OD)^2))
  r2_pool   <- 1 - sum((od_all$OD - od_all$pred_OD)^2) /
    sum((od_all$OD - mean(od_all$OD))^2)

  p_loso <- ggplot(od_all, aes(pred_OD, OD, colour = held_out, shape = Diversity)) +
    geom_abline(linetype = 2) +
    geom_point(alpha = 0.9, size = 2) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("RMSE = %.2f\nR² = %.2f", rmse_pool, r2_pool),
             hjust = -0.05, vjust = 1.4, size = 4) +
    scale_colour_manual(values = pal_short, name = "Held-out stress") +
    scale_shape_manual(values = c("2" = 16, "4" = 17, "8" = 15), name = "Richness") +
    labs(x = expression(Predicted~OD[600]~(held-out~regime)),
         y = expression(OD[600])) +
    theme_classic(base_family = "Helvetica") +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))

  ggsave(P_FIG("Fig_S8_losocv.png"),  p_loso, width = 170, height = 130, units = "mm", dpi = 300)
  ggsave(P_FIG("Fig_S8_losocv.tiff"), p_loso, width = 170, height = 130, units = "mm",
         dpi = 600, device = ragg::agg_tiff, compression = "lzw")
}

print(by_stress)
print(summary_tbl)
message("Done: Tables S20-S21 and Fig. S8 written.")
