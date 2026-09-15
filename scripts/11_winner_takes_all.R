# ======================================================================
# 11_winner_takes_all.R
#
# Winner-takes-all baseline (Reviewer 4, Ecology Letters review 2026):
# how well can composition and biomass be predicted by simply assuming
# that the community member with the highest monoculture growth rate in
# a given stress regime dominates completely?
#
# For each community sample, the predicted composition assigns relative
# abundance 1 to the member with the highest raw monoculture growth rate
# under that regime and 0 to all others. The corresponding "AWM" is then
# just the winner's growth rate, and biomass is predicted with the same
# regression used for the other baselines (OD ~ AWM * Stress + Diversity,
# in-sample, as in Table S11).
#
# No Stan fitting required. Runtime: seconds.
#
# Outputs:
#   results/tables/Table_S25_winner_takes_all.csv
# ======================================================================

source("scripts/utils_bayes_prep.R")

message("[1/3] Preparing data ...")
trait_r <- growth_table_default()
prep    <- prep_bayes_data(trait_r)
comm_od <- load_comm_od()

ord     <- prep$ord_full
SAMPLES <- prep$SAMPLES

# ---------------------- Winner-takes-all composition -----------------
message("[2/3] Building winner-takes-all predictions ...")

wta <- ord %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(p_hat = as.numeric(g == max(g, na.rm = TRUE)),
                p_hat = p_hat / sum(p_hat)) %>%   # split ties evenly
  dplyr::ungroup() %>%
  dplyr::rename(p_obs = Abundance)

# Composition accuracy: weighted metrics (as Fig 3a) + per-sample means
wm <- weighted_comp_metrics(wta %>% dplyr::select(p_obs, p_hat))

comp_per_sample <- wta %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    RMSE = sqrt(mean((p_obs - p_hat)^2)),
    JS   = jsd_local(p_obs, p_hat),
    .groups = "drop"
  )

# How often is the predicted winner actually the observed dominant taxon?
winner_hit <- wta %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    hit = Id[which.max(p_hat)] == Id[which.max(p_obs)],
    .groups = "drop"
  )

# ---------------------- Biomass via winner growth rate ---------------
message("[3/3] Biomass regression ...")

awm_wta <- wta %>%
  dplyr::group_by(SampleID, Com_Id, Stress, Diversity) %>%
  dplyr::summarise(AWM = sum(p_hat * g), .groups = "drop") %>%
  dplyr::group_by(Com_Id, Stress, Diversity) %>%
  dplyr::summarise(AWM = mean(AWM), .groups = "drop") %>%
  dplyr::mutate(Stress = as.character(Stress), Diversity = as.character(Diversity))

df_wta <- comm_od %>%
  dplyr::inner_join(awm_wta, by = c("Com_Id","Stress","Diversity")) %>%
  dplyr::filter(is.finite(OD), is.finite(AWM)) %>%
  dplyr::mutate(Stress = factor(Stress), Diversity = factor(Diversity))

m_wta <- lm(OD ~ AWM * Stress + Diversity, data = df_wta)
pred  <- predict(m_wta, newdata = df_wta)

od_r2   <- 1 - sum((df_wta$OD - pred)^2) / sum((df_wta$OD - mean(df_wta$OD))^2)
od_rmse <- sqrt(mean((df_wta$OD - pred)^2))

tab <- tibble::tibble(
  model            = "Winner-takes-all (highest growth rate dominates)",
  comp_wR2         = round(wm["wR2"], 3),
  comp_wRMSE       = round(wm["wRMSE"], 3),
  comp_RMSE_mean   = round(mean(comp_per_sample$RMSE), 3),
  comp_JS_mean     = round(mean(comp_per_sample$JS), 3),
  winner_correct   = round(mean(winner_hit$hit), 3),
  OD_R2            = round(od_r2, 3),
  OD_RMSE          = round(od_rmse, 3),
  n_samples        = nrow(SAMPLES),
  n_OD             = nrow(df_wta)
)

readr::write_csv(tab, P_TAB("Table_S25_winner_takes_all.csv"))
print(as.data.frame(tab))
message("Done: Table S25 written.")
