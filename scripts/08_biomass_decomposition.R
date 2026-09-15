# ======================================================================
# 08_biomass_decomposition.R
#
# Marginal contribution of the growth trait to the biomass (OD600) model,
# as requested by Reviewer 3 (Ecology Letters review, 2026) and reported
# in the appeal letter. Produces Table S19.
#
# Compares three nested models for community biomass:
#   Null     : OD ~ Stress + Diversity              (10 parameters)
#   Additive : OD ~ AWM + Stress + Diversity        (11 parameters; Eq. 6 additive)
#   Full     : OD ~ AWM * Stress + Diversity        (18 parameters; model as fitted)
#
# AWM is the abundance-weighted mean growth from the main Bayesian fit
# (results/rds/bayes_fit_partA.rds). Point estimates use the posterior-MI
# procedure of script 05 (pooled predictions across 400 posterior draws);
# uncertainty intervals are the 5th-95th percentiles of the per-draw fit
# metrics across the same draws.
#
# Outputs:
#   results/tables/Table_S19_biomass_model_decomposition.csv
#   results/tables/Table_S19b_trait_term_tests.csv
#   results/rds/decomposition_mi_draws.rds
#
# Requires: results/rds/bayes_fit_partA.rds (run script 05 Part A first).
# Runtime: a few minutes (no Stan refits).
# ======================================================================

source("scripts/utils_bayes_prep.R")

NDRAWS_MI <- 400

# ---------------------- Data & cached fit ----------------------------
message("[1/4] Preparing data ...")
trait_r <- growth_table_default()
prep    <- prep_bayes_data(trait_r)

obj <- build_ragged(prep$ord_full, prep$SAMPLES$SampleID,
                    prep$taxa_levels, prep$stress_lvls, prep$rich_lvls,
                    USE_RICHNESS_KAPPA = FALSE, use_taxon_biases_int = 1L)

comm_od <- load_comm_od()

fit_cache <- P_RDS("bayes_fit_partA.rds")
stopifnot(file.exists(fit_cache))
message("[2/4] Loading cached Stan fit ...")
fit   <- readRDS(fit_cache)
draws <- rstan::extract(fit)
nd    <- if (!is.null(draws$kappa)) dim(draws$kappa)[1] else length(draws$mu_kappa)

set.seed(123)
sel_mi <- if (nd > NDRAWS_MI) sample.int(nd, NDRAWS_MI) else seq_len(nd)

# ---------------------- Per-draw AWM ---------------------------------
message("[3/4] Computing per-draw AWM and refitting OD models ...")

samples   <- obj$samples_sub
idx_split <- obj$idx_split
g_z_vec   <- as.numeric(obj$ord_sub$g_z)
g_raw_vec <- as.numeric(obj$ord_sub$g)

awm_one_draw <- function(d) {
  A <- numeric(nrow(samples))
  for (n in seq_len(nrow(samples))) {
    seg <- idx_split[[n]]
    kap <- get_kappa_cv(draws, d, obj$stress_id[n], 1L)
    p   <- softmax_vec(kap * g_z_vec[seg] + draws$delta[d, obj$tax_of[seg]])
    A[n] <- sum(p * g_raw_vec[seg])
  }
  tibble::tibble(SampleID = samples$SampleID, AWM = A) %>%
    dplyr::left_join(samples, by = "SampleID") %>%
    dplyr::group_by(Com_Id, Stress, Diversity) %>%
    dplyr::summarise(AWM = mean(AWM), .groups = "drop") %>%
    dplyr::mutate(Stress = as.character(Stress), Diversity = as.character(Diversity))
}

f_full <- OD ~ AWM * Stress + Diversity
f_add  <- OD ~ AWM + Stress + Diversity
f_null <- OD ~ Stress + Diversity

fit_metrics <- function(m, df) {
  pred <- as.numeric(predict(m, newdata = df))
  c(
    R2   = 1 - sum((df$OD - pred)^2) / sum((df$OD - mean(df$OD))^2),
    RMSE = sqrt(mean((df$OD - pred)^2))
  )
}

draw_stats <- vector("list", length(sel_mi))
pred_full_list <- vector("list", length(sel_mi))
pred_add_list  <- vector("list", length(sel_mi))
awm_draw_list  <- vector("list", length(sel_mi))

for (ii in seq_along(sel_mi)) {
  d     <- sel_mi[ii]
  awm_draw_list[[ii]] <- awm_one_draw(d)
  dat_d <- comm_od %>%
    dplyr::left_join(awm_draw_list[[ii]], by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::filter(is.finite(OD), is.finite(AWM)) %>%
    dplyr::mutate(Stress = factor(Stress), Diversity = factor(Diversity))

  m_full_d <- lm(f_full, data = dat_d)
  m_add_d  <- lm(f_add,  data = dat_d)
  m_null_d <- lm(f_null, data = dat_d)

  mf <- fit_metrics(m_full_d, dat_d)
  ma <- fit_metrics(m_add_d,  dat_d)

  rss_null <- sum(residuals(m_null_d)^2)
  rss_add  <- sum(residuals(m_add_d)^2)
  rss_full <- sum(residuals(m_full_d)^2)
  aov_fn   <- anova(m_null_d, m_full_d)

  draw_stats[[ii]] <- tibble::tibble(
    draw       = d,
    R2_full    = mf["R2"],   RMSE_full = mf["RMSE"],
    R2_add     = ma["R2"],   RMSE_add  = ma["RMSE"],
    partialR2  = 1 - rss_add / rss_null,
    partialR2_full = 1 - rss_full / rss_null,
    F_trait    = aov_fn$F[2],
    df1        = aov_fn$Df[2],
    df2        = aov_fn$Res.Df[2],
    t_AWM_add  = summary(m_add_d)$coefficients["AWM", "t value"]
  )

  key <- dat_d %>% dplyr::select(Com_Id, Stress, Diversity, OD)
  pred_full_list[[ii]] <- key %>% dplyr::mutate(pred = predict(m_full_d, newdata = dat_d))
  pred_add_list[[ii]]  <- key %>% dplyr::mutate(pred = predict(m_add_d,  newdata = dat_d))
}

draw_stats <- dplyr::bind_rows(draw_stats)
saveRDS(draw_stats, P_RDS("decomposition_mi_draws.rds"))

# ---------------------- Pooled-MI point estimates --------------------
pool_mi <- function(lst) {
  dplyr::bind_rows(lst) %>%
    dplyr::group_by(Com_Id, Stress, Diversity, OD) %>%
    dplyr::summarise(pred_mean = mean(pred), .groups = "drop")
}

mi_metrics <- function(pooled) {
  c(
    R2   = 1 - sum((pooled$OD - pooled$pred_mean)^2) / sum((pooled$OD - mean(pooled$OD))^2),
    RMSE = sqrt(mean((pooled$OD - pooled$pred_mean)^2))
  )
}

met_full <- mi_metrics(pool_mi(pred_full_list))
met_add  <- mi_metrics(pool_mi(pred_add_list))

# Null model is draw-invariant: fit once on the point-estimate data frame
dat_pt <- comm_od %>%
  dplyr::left_join(
    dplyr::bind_rows(awm_draw_list) %>%
      dplyr::group_by(Com_Id, Stress, Diversity) %>%
      dplyr::summarise(AWM = mean(AWM), .groups = "drop"),
    by = c("Com_Id","Stress","Diversity")
  ) %>%
  dplyr::filter(is.finite(OD), is.finite(AWM)) %>%
  dplyr::mutate(Stress = factor(Stress), Diversity = factor(Diversity))

m_null <- lm(f_null, data = dat_pt)
m_add  <- lm(f_add,  data = dat_pt)
m_full <- lm(f_full, data = dat_pt)
met_null <- fit_metrics(m_null, dat_pt)

n_par <- function(m) sum(!is.na(coef(m)))

q90 <- function(x) stats::quantile(x, c(0.05, 0.95), na.rm = TRUE)

ci_R2_full   <- q90(draw_stats$R2_full);   ci_RMSE_full <- q90(draw_stats$RMSE_full)
ci_R2_add    <- q90(draw_stats$R2_add);    ci_RMSE_add  <- q90(draw_stats$RMSE_add)

# ---------------------- Table S19 ------------------------------------
message("[4/4] Writing tables ...")

tab_s19 <- tibble::tibble(
  Model      = c("Null (stress + richness only)",
                 "Additive trait term (Eq. 6)",
                 "Full model as fitted"),
  R2         = round(c(met_null["R2"],   met_add["R2"],   met_full["R2"]),   3),
  R2_lwr90   = c(NA, round(ci_R2_add[1], 3),   round(ci_R2_full[1], 3)),
  R2_upr90   = c(NA, round(ci_R2_add[2], 3),   round(ci_R2_full[2], 3)),
  RMSE       = round(c(met_null["RMSE"], met_add["RMSE"], met_full["RMSE"]), 3),
  RMSE_lwr90 = c(NA, round(ci_RMSE_add[1], 3), round(ci_RMSE_full[1], 3)),
  RMSE_upr90 = c(NA, round(ci_RMSE_add[2], 3), round(ci_RMSE_full[2], 3)),
  Parameters = c(n_par(m_null), n_par(m_add), n_par(m_full)),
  n          = nrow(dat_pt)
)

readr::write_csv(tab_s19, P_TAB("Table_S19_biomass_model_decomposition.csv"))

# Point-estimate hypothesis tests (on the MI-mean AWM), with MI intervals
aov_pt <- anova(m_null, m_full)
tab_tests <- tibble::tibble(
  quantity = c("partial_R2_additive_over_null",
               "partial_R2_full_over_null",
               "F_trait_terms_full_vs_null", "df1", "df2", "p_value",
               "t_AWM_additive"),
  point    = c(1 - sum(residuals(m_add)^2)  / sum(residuals(m_null)^2),
               1 - sum(residuals(m_full)^2) / sum(residuals(m_null)^2),
               aov_pt$F[2], aov_pt$Df[2], aov_pt$Res.Df[2],
               aov_pt$`Pr(>F)`[2],
               summary(m_add)$coefficients["AWM", "t value"]),
  lwr90    = c(q90(draw_stats$partialR2)[1], q90(draw_stats$partialR2_full)[1],
               q90(draw_stats$F_trait)[1], NA, NA, NA,
               q90(draw_stats$t_AWM_add)[1]),
  upr90    = c(q90(draw_stats$partialR2)[2], q90(draw_stats$partialR2_full)[2],
               q90(draw_stats$F_trait)[2], NA, NA, NA,
               q90(draw_stats$t_AWM_add)[2])
)

readr::write_csv(tab_tests, P_TAB("Table_S19b_trait_term_tests.csv"))

print(tab_s19)
print(tab_tests)
message("Done: Table S19 written to results/tables/.")
