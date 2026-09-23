# ======================================================================
# 10_trait_sensitivity.R
#
# Two reviewer-requested analyses (Ecology Letters review, 2026), run as
# one suite because they share the monoculture growth-curve refits:
#
# (A) Growth-rate estimation sensitivity (Reviewer 3): re-estimate
#     exponential growth rates with simple log-linear fits over the
#     exponential window of each monoculture OD time series, and compare
#     with the Gompertz estimates used in the pipeline.
#
# (B) Alternative monoculture traits (Reviewer 1): compare the predictive
#     performance of lag time and carrying capacity (log10 nmax), already
#     estimated in the growth-curve fits, against growth rate, by running
#     the full composition -> biomass pipeline with each trait in turn.
#
# Steps:
#   1. Refit Gompertz curves (same model, bounds and nls.multstart
#      settings as scripts/07) to the raw monoculture OD time series in
#      data/Cut_OD_data/, restricted to the 8 factorial regimes, to get
#      per-curve r (mumax), lag and log10_nmax; also fit the log-linear
#      exponential-window slope per curve.
#   2. Compare growth-rate estimates (log-linear vs Gompertz refit vs
#      pipeline values) -> Table S23, Fig. S9.
#   3. For each trait (r pipeline [cached fit], r log-linear, lag,
#      log10_nmax): fit the Dirichlet-softmax composition model, then the
#      posterior-MI OD regression -> Table S22.
#
# Outputs:
#   results/tables/Table_S22_trait_comparison.csv
#   results/tables/Table_S23_growthrate_method_comparison.csv
#   results/figures/Fig_S9_loglinear_vs_gompertz.png / .tiff
#   results/rds/trait_refit_curves.rds
#   results/rds/bayes_fit_trait_<name>.rds  (per-trait Stan fits, cached)
#
# Runtime: growth-curve refits ~10-20 min; 3 Stan fits (3000 iter each).
# ======================================================================

source("scripts/utils_bayes_prep.R")
source("scripts/06_growth_curve_models.R")
for (p in c("nls.multstart","ggplot2","ragg")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ggplot2); library(ragg); library(nls.multstart)
})

CHAINS <- 4; ITER <- 3000; WARMUP <- 1500
ADAPT  <- list(adapt_delta = 0.99, max_treedepth = 15)
CORES  <- 1          # set to 4 for parallel chains if stable on your machine
NDRAWS <- 400
OD_FORMULA <- stats::as.formula("OD ~ AWM * Stress + Diversity")

# ======================================================================
# 1) LOAD RAW MONOCULTURE TIME SERIES, MAP TO THE 8 REGIMES
# ======================================================================
message("[1/5] Loading raw OD time series ...")

od_floor <- function(x) pmax(as.numeric(x), 0.040)   # as in scripts/07

raw_temp <- readr::read_csv(P_IN("Cut_OD_data/cutPoint_D17_edit2_temp.csv"),
                            show_col_types = FALSE) %>%
  dplyr::filter(TempLevel %in% c(20, 38), Id != "BLANK") %>%
  dplyr::mutate(
    Stress = dplyr::if_else(TempLevel == 38, paste0(Stress, "Temp"), Stress),
    Stress = dplyr::case_when(
      Stress == "ControlTemp" ~ "Temp",
      Stress == "pH_SalTemp"  ~ "pHSalTemp",
      Stress == "pH_Sal"      ~ "pHSal",
      TRUE                    ~ Stress
    ),
    source = "temp"
  ) %>%
  dplyr::transmute(source, Id = trimf(Id), Stress = recode_stress(Stress),
                   Rep = as.character(Rep), t = as.numeric(t), OD = as.numeric(OD))

raw_pH <- readr::read_csv(P_IN("Cut_OD_data/cutPoint_D17_edit2_pH.csv"),
                          show_col_types = FALSE) %>%
  dplyr::filter(Id != "BLANK", pH_Level %in% c(5.5, 7.2)) %>%
  dplyr::mutate(
    Stress = dplyr::case_when(
      (SalinityLevel == 0  & TempLevel == 20) ~ "Control",
      (SalinityLevel == 0  & TempLevel == 38) ~ "Temp",
      (SalinityLevel == 20 & TempLevel == 20) ~ "Sal",
      (SalinityLevel == 20 & TempLevel == 38) ~ "SalTemp"
    )
  ) %>%
  dplyr::filter(!is.na(Stress)) %>%
  dplyr::mutate(
    Stress = dplyr::if_else(pH_Level == 5.5, paste0("pH", Stress), Stress),
    Stress = dplyr::if_else(Stress == "pHControl", "pH", Stress),
    source = "pH"
  ) %>%
  dplyr::transmute(source, Id = trimf(Id), Stress = recode_stress(Stress),
                   Rep = as.character(Rep), t = as.numeric(t), OD = as.numeric(OD))

raw_sal <- readr::read_csv(P_IN("Cut_OD_data/cutPoint_D17_edit3_Sal.csv"),
                           show_col_types = FALSE) %>%
  dplyr::filter(Id != "BLANK", SalinityLevel %in% c(0, 20)) %>%
  dplyr::mutate(
    Stress = dplyr::case_when(
      (pH_Level == 7.2 & TempLevel == 20) ~ "Control",
      (pH_Level == 7.2 & TempLevel == 38) ~ "Temp",
      (pH_Level == 5.5 & TempLevel == 20) ~ "pH",
      (pH_Level == 5.5 & TempLevel == 38) ~ "pHTemp"
    )
  ) %>%
  dplyr::filter(!is.na(Stress)) %>%
  dplyr::mutate(
    Stress = dplyr::if_else(SalinityLevel == 20, paste0(Stress, "Sal"), Stress),
    Stress = dplyr::case_when(
      Stress == "ControlSal" ~ "Sal",
      Stress == "pHTempSal"  ~ "pHSalTemp",
      Stress == "TempSal"    ~ "SalTemp",
      TRUE                   ~ Stress
    ),
    source = "sal"
  ) %>%
  dplyr::transmute(source, Id = trimf(Id), Stress = recode_stress(Stress),
                   Rep = as.character(Rep), t = as.numeric(t), OD = as.numeric(OD))

raw_all <- dplyr::bind_rows(raw_temp, raw_pH, raw_sal) %>%
  dplyr::filter(is.finite(t), is.finite(OD)) %>%
  dplyr::mutate(log10_od_cor = log10(od_floor(OD)))

# ======================================================================
# 2) PER-CURVE FITS: GOMPERTZ (as in 07) + LOG-LINEAR WINDOW
# ======================================================================
message("[2/5] Fitting Gompertz + log-linear per curve (this takes a while) ...")

refit_cache <- P_RDS("trait_refit_curves.rds")
if (file.exists(refit_cache)) {
  message("  Loading cached curve refits: ", refit_cache)
  curve_fits <- readRDS(refit_cache)
} else {

  fit_gompertz_curve <- function(df) {
    m <- tryCatch(
      nls_multstart(
        log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
        data = df, iter = 500,
        start_lower = c(log10_nmax = -3,   log10_n0 = -3, mumax = 0,  lag = 0),
        start_upper = c(log10_nmax = -0.2, log10_n0 = -1, mumax = 5,  lag = 1500),
        supp_errors = "Y", na.action = na.omit,
        lower =       c(log10_nmax = -5,   log10_n0 = -5, mumax = 0,  lag = 0),
        upper =       c(log10_nmax =  0,   log10_n0 =  0, mumax = 10, lag = 3000)
      ),
      error = function(e) NULL
    )
    if (is.null(m)) return(c(r_gomp = NA, lag = NA, log10_nmax = NA))
    cf <- coef(m)
    c(r_gomp = unname(cf["mumax"]), lag = unname(cf["lag"]),
      log10_nmax = unname(cf["log10_nmax"]))
  }

  # Log-linear: steepest significant slope of log(OD) over a rolling
  # window of consecutive time points (the exponential window).
  fit_loglinear_curve <- function(df, wmin = 4) {
    df <- df %>% dplyr::arrange(t)
    n  <- nrow(df)
    if (n < wmin + 1) return(NA_real_)
    best <- NA_real_
    for (w in wmin:min(n, 8)) {
      for (i in seq_len(n - w + 1)) {
        sub <- df[i:(i + w - 1), ]
        if (diff(range(sub$t)) <= 0) next
        sl <- tryCatch(coef(lm(log10_od_cor ~ t, data = sub))[2], error = function(e) NA)
        if (is.finite(sl) && (is.na(best) || sl > best)) best <- sl
      }
    }
    # convert log10 slope to ln-based rate, comparable with Gompertz mumax
    unname(best) * log(10)
  }

  curve_fits <- raw_all %>%
    dplyr::group_by(source, Id, Stress, Rep) %>%
    dplyr::group_modify(function(df, key) {
      if (nrow(df) < 8) return(tibble::tibble(
        r_gomp = NA_real_, lag = NA_real_, log10_nmax = NA_real_, r_loglin = NA_real_
      ))
      g  <- fit_gompertz_curve(df)
      ll <- fit_loglinear_curve(df)
      tibble::tibble(
        r_gomp = g["r_gomp"], lag = g["lag"],
        log10_nmax = g["log10_nmax"], r_loglin = ll
      )
    }) %>%
    dplyr::ungroup()

  saveRDS(curve_fits, refit_cache)
  message("  Saved curve refits: ", refit_cache)
}

# ======================================================================
# 3) GROWTH-RATE METHOD COMPARISON (Reviewer 3) -> Table S23, Fig S9
# ======================================================================
message("[3/5] Comparing growth-rate estimation methods ...")

mean_by <- function(df, col) {
  df %>%
    dplyr::group_by(Id, Stress) %>%
    dplyr::summarise(val = mean(.data[[col]], na.rm = TRUE), .groups = "drop")
}

r_pipeline <- load_growth_replicates() %>% mean_by("estimate") %>% dplyr::rename(r_pipeline = val)
r_gomp_ref <- mean_by(curve_fits, "r_gomp")   %>% dplyr::rename(r_gompertz_refit = val)
r_loglin_m <- mean_by(curve_fits, "r_loglin") %>% dplyr::rename(r_loglinear = val)

r_cmp <- r_pipeline %>%
  dplyr::inner_join(r_gomp_ref, by = c("Id","Stress")) %>%
  dplyr::inner_join(r_loglin_m, by = c("Id","Stress"))

cor_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  c(pearson  = cor(x[ok], y[ok]),
    spearman = cor(x[ok], y[ok], method = "spearman"),
    n        = sum(ok))
}

# Within-stress rank agreement (what the softmax model actually uses)
rank_by_stress <- r_cmp %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(
    spearman_loglin_vs_pipeline = cor(r_loglinear, r_pipeline, method = "spearman",
                                      use = "complete.obs"),
    .groups = "drop"
  )

c1 <- cor_stats(r_cmp$r_loglinear, r_cmp$r_pipeline)
c2 <- cor_stats(r_cmp$r_loglinear, r_cmp$r_gompertz_refit)
c3 <- cor_stats(r_cmp$r_gompertz_refit, r_cmp$r_pipeline)

tab_s23 <- tibble::tibble(
  comparison = c("log-linear vs pipeline Gompertz",
                 "log-linear vs Gompertz refit",
                 "Gompertz refit vs pipeline Gompertz"),
  pearson_r  = round(c(c1["pearson"],  c2["pearson"],  c3["pearson"]),  3),
  spearman_rho = round(c(c1["spearman"], c2["spearman"], c3["spearman"]), 3),
  n_Id_x_Stress = c(c1["n"], c2["n"], c3["n"])
)
readr::write_csv(tab_s23, P_TAB("Table_S23_growthrate_method_comparison.csv"))
readr::write_csv(rank_by_stress, P_TAB("Table_S23b_loglin_rank_agreement_by_stress.csv"))

p_s9 <- ggplot(r_cmp, aes(r_pipeline, r_loglinear, colour = Stress)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.9, size = 2) +
  scale_colour_manual(values = pal_short, name = "Stress") +
  labs(x = expression(Gompertz~growth~rate~(pipeline)),
       y = expression(Log-linear~growth~rate)) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pearson r = %.2f\nSpearman ρ = %.2f",
                           c1["pearson"], c1["spearman"]),
           hjust = -0.05, vjust = 1.4, size = 4) +
  theme_classic(base_family = "Helvetica") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))

ggsave(P_FIG("Fig_S9_loglinear_vs_gompertz.png"),  p_s9, width = 170, height = 130,
       units = "mm", dpi = 300)
ggsave(P_FIG("Fig_S9_loglinear_vs_gompertz.tiff"), p_s9, width = 170, height = 130,
       units = "mm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")

# ======================================================================
# 4) TRAIT-SUBSTITUTED PIPELINE RUNS (Reviewers 1 & 3) -> Table S22
# ======================================================================
message("[4/5] Running composition -> biomass pipeline per trait ...")

# Compiled lazily: sm is a memoising function that compiles the model on
# first call; fit_softmax_cached() calls it only when a cached fit is missing.
sm <- local({ m <- NULL; function() { if (is.null(m)) m <<- compile_softmax(P_STAN("softmax_dirichlet_refit.stan")); m } })
comm_od <- load_comm_od()

# For lag, failed 'without-lag' style fits are structurally lag = 0;
# keep NA (dropped in averaging) since we refit Gompertz-with-lag only.
trait_tables <- list(
  r_gompertz  = growth_table_default(),                    # pipeline trait (cached fit)
  r_loglinear = make_trait_table(curve_fits, "r_loglin"),
  lag_time    = make_trait_table(curve_fits, "lag"),
  carrying_capacity_log10nmax = make_trait_table(curve_fits, "log10_nmax")
)

run_pipeline_for_trait <- function(name, trait_tbl) {
  message("  Trait: ", name)
  prep <- prep_bayes_data(trait_tbl)
  obj  <- build_ragged(prep$ord_full, prep$SAMPLES$SampleID, prep$taxa_levels,
                       prep$stress_lvls, prep$rich_lvls, FALSE, 1L)

  if (name == "r_gompertz" && file.exists(P_RDS("bayes_fit_partA.rds"))) {
    fit <- readRDS(P_RDS("bayes_fit_partA.rds"))
  } else {
    fit <- fit_softmax_cached(
      sm, obj$stan_data,
      cache_path = P_RDS(sprintf("bayes_fit_trait_%s.rds", name)),
      chains = CHAINS, iter = ITER, warmup = WARMUP,
      adapt = ADAPT, seed = 123, cores = CORES
    )
  }

  draws <- rstan::extract(fit)
  nd    <- dim(draws$kappa)[1]
  set.seed(123)
  sel <- if (nd > NDRAWS) sample.int(nd, NDRAWS) else seq_len(nd)

  pred <- predict_comp_and_awm(obj, draws, sel, FALSE, prep$taxa_levels)
  wm   <- weighted_comp_metrics(pred$pred_comp)

  # Posterior-MI OD regression, pooled predictions (as in script 05)
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

  preds <- vector("list", length(sel))
  for (ii in seq_along(sel)) {
    dat_d <- comm_od %>%
      dplyr::inner_join(awm_one_draw(sel[ii]), by = c("Com_Id","Stress","Diversity")) %>%
      dplyr::filter(is.finite(OD), is.finite(AWM)) %>%
      dplyr::mutate(Stress = factor(Stress), Diversity = factor(Diversity))
    if (nrow(dat_d) < 5) next
    m_d <- tryCatch(lm(OD_FORMULA, data = dat_d), error = function(e) NULL)
    if (is.null(m_d)) next
    preds[[ii]] <- dat_d %>%
      dplyr::transmute(Com_Id, Stress, Diversity, OD,
                       pred = as.numeric(predict(m_d, newdata = dat_d)))
  }

  pooled <- dplyr::bind_rows(preds) %>%
    dplyr::group_by(Com_Id, Stress, Diversity, OD) %>%
    dplyr::summarise(pred_mean = mean(pred), .groups = "drop")

  od_r2   <- 1 - sum((pooled$OD - pooled$pred_mean)^2) /
    sum((pooled$OD - mean(pooled$OD))^2)
  od_rmse <- sqrt(mean((pooled$OD - pooled$pred_mean)^2))

  tibble::tibble(
    trait     = name,
    n_samples = nrow(samples),
    comp_wR2  = round(wm["wR2"], 3),
    comp_wRMSE = round(wm["wRMSE"], 3),
    OD_R2     = round(od_r2, 3),
    OD_RMSE   = round(od_rmse, 3)
  )
}

tab_s22 <- dplyr::bind_rows(
  purrr::imap(trait_tables, function(tbl, nm) run_pipeline_for_trait(nm, tbl))
)

readr::write_csv(tab_s22, P_TAB("Table_S22_trait_comparison.csv"))

message("[5/5] Done.")
print(tab_s23)
print(tab_s22)
message("Tables S22-S23 and Fig. S9 written.")
