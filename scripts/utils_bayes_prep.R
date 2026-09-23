# ======================================================================
# utils_bayes_prep.R
#
# Shared data preparation and helpers for the reviewer-response analyses
# (scripts 08-10). Mirrors, line for line where possible, the data prep
# in scripts/05_fig3_bayes_od.R (Parts A and B) so that results are
# directly comparable with the published pipeline.
#
# Exposes:
# - design_communities()          planned community membership
# - load_abundance()              endpoint compositions (zeros inserted)
# - load_growth_replicates()      per-replicate Gompertz growth rates,
#                                 mapped to the 8 factorial stress regimes
# - make_trait_table(df, col)     mean trait per Id x Stress + within-stress z
# - load_comm_od()                community OD600 (t = 2, no Evo)
# - prep_bayes_data(trait_tbl)    aligned ragged structures for Stan
# - build_ragged(), softmax_vec(), get_kappa_cv(), predict_comp_and_awm()
# - fit_softmax_cached()          compile once, fit with RDS caching
# - stan_softmax_code()           the Dirichlet-softmax model code (as in 05)
# ======================================================================

source("scripts/utils_functions.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(readr); library(purrr); library(rstan)
})

# ---------------------- DESIGN (planned members) ---------------------
design_communities <- function() {
  tibble::tribble(
    ~Com_Id, ~Id,
    "Com1","D14", "Com1","I15",
    "Com2","D17", "Com2","I9",
    "Com3","I8",  "Com3","I11",
    "Com4","I2",  "Com4","I18",
    "Com5","I20", "Com5","I23",
    "Com6","I22", "Com6","D11",
    # 4-taxon sets
    "Com7","I20","Com7","I22","Com7","I23","Com7","D14",
    "Com8","I2","Com8","I9","Com8","I15","Com8","D17",
    "Com9","I8","Com9","I11","Com9","I18","Com9","D14",
    "Com10","I2","Com10","I8","Com10","I15","Com10","I23",
    "Com11","I9","Com11","I11","Com11","I20","Com11","D11",
    "Com12","I18","Com12","I22","Com12","D11","Com12","D17",
    # 8-taxon sets
    "Com13","I2","Com13","I9","Com13","I11","Com13","I20","Com13","I22","Com13","D11","Com13","D14","Com13","D17",
    "Com14","I2","Com14","I8","Com14","I9","Com14","I18","Com14","I20","Com14","I22","Com14","D11","Com14","D17",
    "Com15","I9","Com15","I11","Com15","I15","Com15","I18","Com15","I20","Com15","I23","Com15","D11","Com15","D14",
    "Com16","I11","Com16","I18","Com16","I20","Com16","I22","Com16","I23","Com16","D11","Com16","D14","Com16","D17",
    "Com17","I2","Com17","I8","Com17","I9","Com17","I15","Com17","I20","Com17","I23","Com17","D11","Com17","D17",
    "Com18","I2","Com18","I9","Com18","I15","Com18","I18","Com18","I20","Com18","I22","Com18","D11","Com18","D17",
    "Com19","I2","Com19","I8","Com19","I9","Com19","I15","Com19","I18","Com19","I20","Com19","I23","Com19","D14",
    "Com20","I2","Com20","I8","Com20","I9","Com20","I11","Com20","I18","Com20","I20","Com20","I22","Com20","D14",
    "Com21","I2","Com21","I8","Com21","I11","Com21","I15","Com21","I18","Com21","I20","Com21","I23","Com21","D14",
    "Com22","I8","Com22","I11","Com22","I15","Com22","I20","Com22","I22","Com22","I23","Com22","D11","Com22","D17"
  ) %>%
    dplyr::mutate(
      Com_Id = gsub("\\s+","", Com_Id),
      Id     = trimf(Id)
    )
}

recode_stress <- function(x) {
  dplyr::case_when(
    x == "pH_Sal"       ~ "pHSal",
    x == "pH_Sal_Temp"  ~ "pHSalTemp",
    x == "pH_Temp"      ~ "pHTemp",
    x == "Sal_Temp"     ~ "SalTemp",
    x == "Temp_Sal"     ~ "SalTemp",
    TRUE                ~ x
  )
}

# Resolve a data file case-insensitively (the abundance file is named
# "ALl_Abundance_Data.csv" on disk; macOS filesystems are usually
# case-insensitive but Linux is not).
resolve_input <- function(fname) {
  p <- P_IN(fname)
  if (file.exists(p)) return(p)
  hits <- list.files(DIR_INPUT, full.names = TRUE)
  hit  <- hits[tolower(basename(hits)) == tolower(fname)]
  if (length(hit) == 1) return(hit)
  stop("Input file not found: ", fname)
}

# ---------------------- ABUNDANCE ------------------------------------
load_abundance <- function(design = design_communities()) {
  abund_long <- utils::read.csv(resolve_input("All_Abundance_Data.csv")) %>%
    dplyr::transmute(
      SampleID  = trimf(SampleID),
      Com_Id    = gsub("\\s+", "", trimf(Com_Id)),
      Stress    = trimf(Stress),
      Diversity = trimf(Diversity),
      Id        = trimf(Id),
      Abundance = as.numeric(Abundance)
    ) %>%
    dplyr::mutate(
      Diversity = dplyr::case_when(
        Diversity %in% c("One","one")     ~ "1",
        Diversity %in% c("Two","two")     ~ "2",
        Diversity %in% c("Four","four")   ~ "4",
        Diversity %in% c("Eight","eight") ~ "8",
        TRUE                              ~ as.character(Diversity)
      ),
      Diversity = as.numeric(Diversity),
      Stress    = recode_stress(Stress)
    ) %>%
    dplyr::filter(Diversity %in% c(2,4,8)) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Abundance = Abundance / sum(Abundance)) %>%
    dplyr::ungroup()

  abund_long %>%
    dplyr::distinct(SampleID, Com_Id, Stress, Diversity) %>%
    dplyr::left_join(design, by = "Com_Id", relationship = "many-to-many") %>%
    dplyr::left_join(
      abund_long,
      by = c("SampleID","Com_Id","Stress","Diversity","Id"),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(Abundance = tidyr::replace_na(Abundance, 0)) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Abundance = Abundance / sum(Abundance)) %>%
    dplyr::ungroup()
}

# ---------------------- GROWTH (per replicate) -----------------------
# Returns one row per Id x Stress x replicate with the Gompertz growth
# rate ('estimate'), mapped onto the 8 factorial regimes exactly as in 05.
load_growth_replicates <- function() {

  load_growth_file <- function(path, mutate_steps) mutate_steps(utils::read.csv(path))

  mut_temp <- function(d) {
    d %>%
      dplyr::filter(TempLevel %in% c(20, 38)) %>%
      dplyr::mutate(Stress = dplyr::if_else(TempLevel == 38, paste0(Stress, "Temp"), Stress)) %>%
      dplyr::mutate(
        Stress = dplyr::case_when(
          Stress == "ControlTemp" ~ "Temp",
          Stress == "pH_SalTemp"  ~ "pHSalTemp",
          Stress == "pH_Sal"      ~ "pHSal",
          TRUE                    ~ Stress
        )
      ) %>%
      dplyr::select(Id, estimate, Stress)
  }

  mut_pH <- function(d) {
    d %>%
      dplyr::mutate(
        Stress = dplyr::case_when(
          (SalinityLevel == 0  & TempLevel == 20) ~ "Control",
          (SalinityLevel == 0  & TempLevel == 38) ~ "Temp",
          (SalinityLevel == 20 & TempLevel == 20) ~ "Sal",
          (SalinityLevel == 20 & TempLevel == 38) ~ "SalTemp"
        )
      ) %>%
      dplyr::filter(!is.na(Stress)) %>%
      dplyr::filter(pH_Level %in% c(5.5, 7.2)) %>%
      dplyr::mutate(Stress = dplyr::if_else(pH_Level == 5.5, paste0("pH", Stress), Stress)) %>%
      dplyr::mutate(Stress = dplyr::if_else(Stress == "pHControl", "pH", Stress)) %>%
      dplyr::select(Id, estimate, Stress)
  }

  mut_sal <- function(d) {
    d %>%
      dplyr::mutate(
        Stress = dplyr::case_when(
          (pH_Level == 7.2 & TempLevel == 20) ~ "Control",
          (pH_Level == 7.2 & TempLevel == 38) ~ "Temp",
          (pH_Level == 5.5 & TempLevel == 20) ~ "pH",
          (pH_Level == 5.5 & TempLevel == 38) ~ "pHTemp"
        )
      ) %>%
      dplyr::filter(!is.na(Stress)) %>%
      dplyr::filter(SalinityLevel %in% c(0, 20)) %>%
      dplyr::mutate(Stress = dplyr::if_else(SalinityLevel == 20, paste0(Stress, "Sal"), Stress)) %>%
      dplyr::mutate(
        Stress = dplyr::case_when(
          Stress == "ControlSal"  ~ "Sal",
          Stress == "pHTempSal"   ~ "pHSalTemp",
          Stress == "TempSal"     ~ "SalTemp",
          TRUE                    ~ Stress
        )
      ) %>%
      dplyr::select(Id, estimate, Stress)
  }

  dplyr::bind_rows(
    load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv"),      mut_temp),
    load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv"), mut_pH),
    load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv"),       mut_sal)
  ) %>%
    dplyr::mutate(Stress = recode_stress(Stress), Id = trimf(Id))
}

# ---------------------- TRAIT TABLE ----------------------------------
# df must have columns Id, Stress and the trait column; replicate values
# are averaged per Id x Stress and z-scored within each stress regime,
# exactly as growth rate is treated in 05.
make_trait_table <- function(df, trait_col = "estimate") {
  df %>%
    dplyr::group_by(Id, Stress) %>%
    dplyr::summarise(g = mean(.data[[trait_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Id = trimf(Id)) %>%
    dplyr::group_by(Stress) %>%
    dplyr::mutate(g_z = as.numeric(scale(g))) %>%
    dplyr::ungroup()
}

# Default: mean Gompertz growth rate per Id x Stress (as in Part A/B)
growth_table_default <- function() make_trait_table(load_growth_replicates())

# ---------------------- COMMUNITY OD ---------------------------------
load_comm_od <- function() {
  utils::read.csv(P_IN("CommunityOD_All_Jan23_edit_zeros.csv")) %>%
    dplyr::mutate(
      Stress    = recode_stress(Stress),
      Com_Id    = paste0("Com", gsub("\\s+","", Community)),
      Diversity = as.character(Diversity_Level)
    ) %>%
    dplyr::filter(Evo_Treatment != "Evo", t == 2) %>%
    dplyr::transmute(
      Com_Id    = gsub("\\s+","", Com_Id),
      Stress    = as.character(Stress),
      Diversity = as.character(Diversity),
      OD        = as.numeric(OD)
    ) %>%
    dplyr::filter(is.finite(OD))
}

# ---------------------- ALIGN & STRUCTURE ----------------------------
# trait_tbl: Id, Stress, g, g_z (from make_trait_table)
prep_bayes_data <- function(trait_tbl, design = design_communities()) {

  abund_complete <- load_abundance(design)

  abund_growth <- abund_complete %>%
    dplyr::left_join(trait_tbl, by = c("Id","Stress"))

  complete_ids <- abund_growth %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(all_g = all(!is.na(g)), .groups = "drop") %>%
    dplyr::filter(all_g) %>%
    dplyr::pull(SampleID)

  ord_full <- abund_growth %>%
    dplyr::filter(SampleID %in% complete_ids) %>%
    dplyr::arrange(SampleID, Id)

  SAMPLES <- ord_full %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      Com_Id    = dplyr::first(Com_Id),
      Stress    = dplyr::first(Stress),
      Diversity = dplyr::first(Diversity),
      len       = dplyr::n(),
      .groups   = "drop"
    )

  taxa_levels <- ord_full %>% dplyr::distinct(Id) %>% dplyr::arrange(Id) %>% dplyr::pull(Id)

  list(
    ord_full    = ord_full,
    SAMPLES     = SAMPLES,
    taxa_levels = taxa_levels,
    stress_lvls = sort(unique(SAMPLES$Stress)),
    rich_lvls   = sort(unique(SAMPLES$Diversity))
  )
}

# ---------------------- RAGGED STRUCTURES (as in 05 Part B) ----------
build_ragged <- function(ord_df, sample_ids_keep,
                         taxa_levels, stress_lvls, rich_lvls,
                         USE_RICHNESS_KAPPA = FALSE, use_taxon_biases_int = 1L) {

  ord_sub <- ord_df %>%
    dplyr::filter(SampleID %in% sample_ids_keep) %>%
    dplyr::arrange(SampleID, Id)

  samples_sub <- ord_sub %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      Com_Id    = dplyr::first(Com_Id),
      Stress    = dplyr::first(Stress),
      Diversity = dplyr::first(Diversity),
      len       = dplyr::n(),
      .groups   = "drop"
    )

  N <- nrow(samples_sub)
  len <- samples_sub$len
  start_idx <- c(1L, 1L + head(cumsum(len), -1L))
  J <- sum(len)

  map_tax <- setNames(seq_along(taxa_levels), taxa_levels)
  tax_of <- as.integer(map_tax[ord_sub$Id])

  map_stress <- setNames(seq_along(stress_lvls), stress_lvls)
  stress_id <- as.integer(map_stress[samples_sub$Stress])

  if (USE_RICHNESS_KAPPA) {
    rich_id <- as.integer(factor(samples_sub$Diversity, levels = rich_lvls))
    Rdim <- length(rich_lvls)
  } else {
    rich_id <- rep(1L, N)
    Rdim <- 1L
  }

  idx_split <- split(seq_len(J), rep.int(seq_len(N), times = len))

  list(
    ord_sub = ord_sub,
    samples_sub = samples_sub,
    idx_split = idx_split,
    tax_of = tax_of,
    stress_id = stress_id,
    rich_id = rich_id,
    stan_data = list(
      N = N, K = length(taxa_levels), J = J, S = length(stress_lvls),
      start_idx = start_idx, len = len,
      p_obs = as.vector(ord_sub$Abundance),
      g_z   = as.vector(ord_sub$g_z),
      tax_of = tax_of,
      stress_id = stress_id,
      Rdim = Rdim,
      rich_id = rich_id,
      use_taxon_biases = use_taxon_biases_int
    )
  )
}

# Jensen-Shannon divergence (log2 units, matching philentropy::JSD default)
jsd_local <- function(p, q) {
  m <- 0.5 * (p + q)
  kl <- function(a, b) { ok <- a > 0; sum(a[ok] * log2(a[ok] / b[ok])) }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

softmax_vec <- function(x) {
  ex <- exp(x - max(x))
  ex / sum(ex)
}

get_kappa_cv <- function(draws, d, s, r = 1L) {
  kd <- draws$kappa
  if (length(dim(kd)) == 2L) kd[d, s] else kd[d, s, r]
}

predict_comp_and_awm <- function(obj, draws, sel_draws, USE_RICHNESS_KAPPA, taxa_levels) {

  idx_split <- obj$idx_split
  tax_of    <- obj$tax_of
  stress_id <- obj$stress_id
  rich_id   <- obj$rich_id
  ord_sub   <- obj$ord_sub
  samples   <- obj$samples_sub

  g_z_vec   <- as.numeric(ord_sub$g_z)
  g_raw_vec <- as.numeric(ord_sub$g)

  pred_rows <- vector("list", nrow(samples))
  awm_rows  <- vector("list", nrow(samples))

  for (n in seq_len(nrow(samples))) {
    seg <- idx_split[[n]]
    s   <- stress_id[n]
    r_i <- if (USE_RICHNESS_KAPPA) rich_id[n] else 1L

    gsub_z  <- g_z_vec[seg]
    gsub_r  <- g_raw_vec[seg]
    taxa_ix <- tax_of[seg]

    P <- matrix(NA_real_, nrow = length(sel_draws), ncol = length(seg))
    A <- numeric(length(sel_draws))

    for (ii in seq_along(sel_draws)) {
      d <- sel_draws[ii]
      kap <- get_kappa_cv(draws, d, s, r_i)
      logits <- kap * gsub_z + draws$delta[d, taxa_ix]
      p <- softmax_vec(logits)
      P[ii, ] <- p
      A[ii]   <- sum(p * gsub_r)
    }

    p_hat <- colMeans(P)

    pred_rows[[n]] <- tibble::tibble(
      SampleID  = samples$SampleID[n],
      Com_Id    = samples$Com_Id[n],
      Stress    = samples$Stress[n],
      Diversity = samples$Diversity[n],
      Id        = taxa_levels[taxa_ix],
      p_hat     = p_hat,
      p_obs     = ord_sub$Abundance[seg]
    )

    awm_rows[[n]] <- tibble::tibble(
      SampleID  = samples$SampleID[n],
      Com_Id    = samples$Com_Id[n],
      Stress    = samples$Stress[n],
      Diversity = samples$Diversity[n],
      AWM_mean  = mean(A, na.rm = TRUE),
      AWM_sd    = stats::sd(A, na.rm = TRUE)
    )
  }

  list(
    pred_comp = dplyr::bind_rows(pred_rows),
    awm_post  = dplyr::bind_rows(awm_rows)
  )
}

# ---------------------- WEIGHTED COMPOSITION METRICS -----------------
# Abundance-weighted RMSE / R2, as in Fig 3a of script 05.
weighted_comp_metrics <- function(pred_comp) {
  w    <- pred_comp$p_obs
  y    <- pred_comp$p_obs
  yhat <- pred_comp$p_hat
  wSSE   <- sum(w * (y - yhat)^2, na.rm = TRUE)
  ybar_w <- sum(w * y, na.rm = TRUE) / sum(w, na.rm = TRUE)
  wSST   <- sum(w * (y - ybar_w)^2, na.rm = TRUE)
  c(wRMSE = sqrt(wSSE / sum(w, na.rm = TRUE)), wR2 = 1 - (wSSE / wSST))
}

# ---------------------- STAN MODEL -----------------------------------
stan_softmax_code <- function() {
'
functions {
  vector segment_vector(vector x, int start, int L) {
    vector[L] y;
    for (m in 1:L) y[m] = x[start + m - 1];
    return y;
  }
}
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> J;
  int<lower=1> S;
  array[N] int<lower=1> start_idx;
  array[N] int<lower=1> len;
  vector[J] p_obs;
  vector[J] g_z;
  array[J] int<lower=1, upper=K> tax_of;
  array[N] int<lower=1, upper=S> stress_id;
  int<lower=1> Rdim;
  array[N] int<lower=1, upper=Rdim> rich_id;
  int<lower=0, upper=1> use_taxon_biases;
}
parameters {
  real mu_kappa;
  real<lower=0> sigma_kappa;
  matrix[S, Rdim] kappa_raw;
  vector[K-1] delta_raw;
  real<lower=0> sigma_delta;
  real mu_phi;
  real<lower=0> sigma_phi;
  vector[S] log_phi_raw;
}
transformed parameters {
  matrix[S, Rdim] kappa;
  vector[K] delta;
  vector[S] log_phi;
  vector[S] phi;

  kappa = mu_kappa + sigma_kappa * kappa_raw;

  delta[1:(K-1)] = use_taxon_biases * sigma_delta * delta_raw;
  delta[K]       = 0;

  log_phi = mu_phi + sigma_phi * log_phi_raw;
  for (s in 1:S) phi[s] = exp(log_phi[s]);
}
model {
  mu_kappa    ~ normal(0, 1);
  sigma_kappa ~ normal(0, 0.5);
  to_vector(kappa_raw) ~ normal(0, 1);

  sigma_delta ~ normal(0, 1);
  delta_raw   ~ normal(0, 1);

  mu_phi    ~ normal(log(50), 1);
  sigma_phi ~ normal(0, 1);
  log_phi_raw ~ normal(0, 1);

  for (n in 1:N) {
    int a = start_idx[n];
    int L = len[n];
    vector[L] eta;
    vector[L] pseg = segment_vector(p_obs, a, L);

    pseg = pseg + 1e-12;
    pseg = pseg / sum(pseg);

    for (m in 1:L) {
      int j = a + m - 1;
      int t = tax_of[j];
      eta[m] = kappa[stress_id[n], rich_id[n]] * g_z[j] + delta[t];
    }

    target += dirichlet_lpdf(pseg | phi[stress_id[n]] * softmax(eta));
  }
}
'
}

# Compile the CV/refit Stan model once per session
compile_softmax <- function(stan_file = P_STAN("softmax_dirichlet_refit.stan")) {
  writeLines(paste0(stan_softmax_code(), "\n"), stan_file)
  rstan::stan_model(stan_file)
}

# Fit with RDS caching (delete the cache file to force a refit)
fit_softmax_cached <- function(sm, stan_data, cache_path,
                               chains = 4, iter = 2000, warmup = 1000,
                               adapt = list(adapt_delta = 0.99, max_treedepth = 15),
                               seed = 123, cores = 1) {
  if (file.exists(cache_path)) {
    message("  Loading cached fit: ", cache_path)
    return(readRDS(cache_path))
  }
  fit <- rstan::sampling(
    sm, data = stan_data,
    chains = chains, iter = iter, warmup = warmup,
    control = adapt, seed = seed, cores = cores
  )
  saveRDS(fit, cache_path)
  message("  Saved fit: ", cache_path)
  fit
}
