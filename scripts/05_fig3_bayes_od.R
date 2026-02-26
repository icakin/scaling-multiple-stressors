# ======================================================================
# Fig 3 + Blocked CV (COMBINED SINGLE SCRIPT)
#
# Part A: Fig 3: Bayesian Composition Model and OD Predictions (main model)
# Part B: Blocked cross-validation (by Com_Id or Stress) + OD fold diagnostic
# ======================================================================

# ------------------------------ PART A -------------------------------
# ======================================================================
# Fig 3: Bayesian Composition Model and OD Predictions
#
# Main figures:
# - Fig 3a: Composition fit (predicted vs observed relative abundances; pooled)
# - Fig 3b: OD predictions using Bayesian AWM × g (posterior MI)
#
# Supplementary figures:
# - Fig S4: Facetted composition fits by Stress (predicted vs observed)
# - Fig S5a: OD predictions using observed AWM × g (oracle ceiling)
# - Fig S5b: OD predictions using equal-abundance × g baseline
# - Fig S5c: OD predictions using trait-shuffled × g baseline
#
# Tables (Part A):
# - Table S9:  Abundance fit metrics (RMSE, R²) by Stress
# - Table S10: elpd_loo by Stress (if log_lik is available)
# - Table S11: Overall OD performance (RMSE, R²) for Bayes + baselines
# - Table S12: JS divergence by Stress
# - Table S13: LOO summary statistics for the Bayes model
#
# Diagnostics (Part A; currently written as "SXX" placeholders):
# - Table S14: Stan key-parameter diagnostics (n_eff, Rhat, etc.)
# - Table SXX: Stan sampler diagnostics (divergences, treedepth hits)
# - Table SXX: Composition 90% interval coverage by Stress
#
# Overview:
# - Loads community composition and single-taxon growth data
# - Fits a Dirichlet–softmax Bayesian model linking z-scored growth (g_z)
#   and (optional) taxon-specific bias terms to observed relative abundances
# - Generates posterior predicted compositions and abundance-weighted mean
#   growth (AWM) per community
# - Predicts OD600 from AWM × g using:
#     * Bayesian AWM (posterior MI)  -> Fig 3b
#     * Observed AWM (oracle)        -> Fig S5a
#     * Equal-abundance baseline     -> Fig S5b
#     * Trait-shuffled baseline      -> Fig S5c
# - Summarises composition and OD prediction performance into Tables S9–S13
#
# Folder conventions:
# - Inputs (composition, growth fits, OD):  data/              (via P_IN())
# - Scripts & helpers:                      scripts/
# - Plots:                                  results/figures/   (via P_FIG())
# - Tables (CSV):                           results/tables/    (via P_TAB())
# - RDS objects (Bayes outputs, AWM, etc.): results/rds/        (via P_RDS())
# ======================================================================

# --- Setup ---
source("scripts/utils_functions.R")
ensure_packages()  # <- no arguments (matches utils_functions.R)

suppressPackageStartupMessages({
  # utils already loads most of these, but calling again is harmless
  library(ggplot2); library(grid)
  library(dplyr);   library(tidyr);  library(tibble)
  library(rstan);   library(loo);    library(philentropy)
  library(ragg);    library(readr); library(purrr)
})

rstan::rstan_options(auto_write = TRUE)
options(stringsAsFactors = FALSE)
set.seed(123)

# ---------------------- USER OPTIONS ---------------------------------
USE_RICHNESS_KAPPA   <- FALSE
USE_TAXON_BIASES     <- TRUE
USE_NESTED_CV_FOR_OD <- FALSE
K_OUTER <- 5; K_INNER <- 3

CHAINS <- 4; ITER <- 3000; WARMUP <- 1500
ADAPT  <- list(adapt_delta = 0.99, max_treedepth = 15)
NDRAWS_COMPOSITION <- 400
NDRAWS_MI          <- 400

# ---------------------- DESIGN (planned members) ---------------------
design <- tibble::tribble(
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

# ---------------------- LOAD & PREP DATA -----------------------------
message("[1/6] Loading & preparing abundance data ...")

abund_long <- read.csv(P_IN("All_Abundance_Data.csv")) %>%
  dplyr::transmute(
    SampleID  = trimf(SampleID),
    Com_Id    = gsub("\\s+", "", trimf(Com_Id)),
    Stress    = trimf(Stress),
    Diversity = trimf(Diversity),
    Id        = trimf(Id),
    Abundance = as.numeric(Abundance)
  ) %>%
  dplyr::mutate(
    Diversity = dplyr::recode(
      Diversity,
      "One"="1","one"="1","Two"="2","two"="2",
      "Four"="4","four"="4","Eight"="8","eight"="8",
      .default = Diversity
    ),
    Diversity = as.numeric(Diversity),
    Stress = dplyr::case_when(
      Stress == "pH_Sal"       ~ "pHSal",
      Stress == "pH_Sal_Temp"  ~ "pHSalTemp",
      Stress == "pH_Temp"      ~ "pHTemp",
      Stress == "Sal_Temp"     ~ "SalTemp",
      TRUE                     ~ Stress
    )
  ) %>%
  dplyr::filter(Diversity %in% c(2,4,8)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Abundance = Abundance / sum(Abundance)) %>%
  dplyr::ungroup()

abund_complete <- abund_long %>%
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

message("[2/6] Loading growth rate files ...")

load_growth_file <- function(path, mutate_steps) {
  mutate_steps(read.csv(path))
}

mut_temp <- function(d) {
  d %>%
    dplyr::filter(TempLevel %in% c(20, 38)) %>%
    dplyr::mutate(Stress = if_else(TempLevel == 38, paste0(Stress, "Temp"), Stress)) %>%
    dplyr::mutate(
      Stress = dplyr::case_when(
        Stress == "ControlTemp" ~ "Temp",
        Stress == "pH_SalTemp"  ~ "pHSalTemp",
        Stress == "pH_Sal"      ~ "pHSal",
        TRUE                    ~ Stress
      )
    ) %>%
    dplyr::select(-TempLevel, -Rep)
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
    dplyr::filter(!is.na(Stress)) %>%   # <-- ADDED: avoid NA Stress rows
    dplyr::select(-SalinityLevel, -TempLevel, -X, -X.1) %>%
    dplyr::filter(pH_Level %in% c(5.5, 7.2)) %>%
    dplyr::mutate(Stress = if_else(pH_Level == 5.5, paste0("pH", Stress), Stress)) %>%
    dplyr::mutate(Stress = if_else(Stress == "pHControl", "pH", Stress)) %>%
    dplyr::select(-pH_Level, -Rep)
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
    dplyr::filter(!is.na(Stress)) %>%   # <-- ADDED: avoid NA Stress rows
    dplyr::select(-pH_Level, -TempLevel, -X) %>%
    dplyr::filter(SalinityLevel %in% c(0, 20)) %>%
    dplyr::mutate(Stress = if_else(SalinityLevel == 20, paste0(Stress, "Sal"), Stress)) %>%
    dplyr::mutate(
      Stress = dplyr::case_when(
        Stress == "ControlSal"  ~ "Sal",
        Stress == "pHTempSal"   ~ "pHSalTemp",
        Stress == "TempSal"     ~ "SalTemp",
        TRUE                    ~ Stress
      )
    ) %>%
    dplyr::select(-SalinityLevel, -Rep)
}

all_growth_raw <- dplyr::bind_rows(
  load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv"), mut_temp),
  load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv"),   mut_pH),
  load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv"),         mut_sal)
)

all_growth <- all_growth_raw %>%
  dplyr::mutate(
    Stress = dplyr::case_when(
      Stress == "pH_Sal"       ~ "pHSal",
      Stress == "pH_Sal_Temp"  ~ "pHSalTemp",
      Stress == "pH_Temp"      ~ "pHTemp",
      Stress == "Temp_Sal"     ~ "SalTemp",
      TRUE                     ~ Stress
    )
  ) %>%
  dplyr::group_by(Id, Stress) %>%
  dplyr::summarise(g = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Id = trimf(Id)) %>%
  dplyr::group_by(Stress) %>%
  dplyr::mutate(g_z = as.numeric(scale(g))) %>%
  dplyr::ungroup()

message("[3/6] Aligning abundance & growth; building ragged structures ...")

abund_growth <- abund_complete %>%
  dplyr::left_join(all_growth, by = c("Id","Stress"))

complete_ids <- abund_growth %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(all_g = all(!is.na(g)), .groups = "drop") %>%
  dplyr::filter(all_g) %>%
  dplyr::pull(SampleID)

data_all <- abund_growth %>%
  dplyr::filter(SampleID %in% complete_ids)

ord <- data_all %>%
  dplyr::arrange(SampleID, Id)

SAMPLES <- ord %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    Com_Id    = dplyr::first(Com_Id),
    Stress    = dplyr::first(Stress),
    Diversity = dplyr::first(Diversity),
    len       = dplyr::n(),
    .groups   = "drop"
  )

N   <- nrow(SAMPLES)
len <- SAMPLES$len
start_idx <- c(1L, 1L + head(cumsum(len), -1L))
J   <- sum(len)

p_obs <- ord$Abundance
g_z   <- ord$g_z
g_raw <- ord$g

K <- ord %>%
  dplyr::distinct(Id) %>%
  dplyr::arrange(Id) %>%
  nrow()

taxa_levels <- ord %>%
  dplyr::distinct(Id) %>%
  dplyr::arrange(Id) %>%
  dplyr::pull(Id)

map_tax <- setNames(seq_along(taxa_levels), taxa_levels)

tax_of      <- as.integer(map_tax[ord$Id])
stress_lvls <- sort(unique(ord$Stress))
S           <- length(stress_lvls)
map_stress  <- setNames(seq_along(stress_lvls), stress_lvls)
stress_id   <- as.integer(map_stress[SAMPLES$Stress])

rich_lvls <- sort(unique(SAMPLES$Diversity))
Rdim      <- if (USE_RICHNESS_KAPPA) length(rich_lvls) else 1L
rich_id   <- if (USE_RICHNESS_KAPPA) {
  as.integer(factor(SAMPLES$Diversity, levels = rich_lvls))
} else {
  rep(1L, N)
}

# sanity: each segment sums to 1
stopifnot(
  all(
    abs(tapply(p_obs, rep(seq_len(N), times = len), sum) - 1) < 1e-6
  )
)

# ---------------------- STAN MODEL -----------------------------------
message("[4/6] Writing, compiling, and sampling Stan model ...")

stan_dirichlet <- '
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
  int<lower=1> start_idx[N];
  int<lower=1> len[N];
  vector[J] p_obs;
  vector[J] g_z;
  int<lower=1, upper=K> tax_of[J];
  int<lower=1, upper=S> stress_id[N];
  int<lower=1> Rdim;
  int<lower=1, upper=Rdim> rich_id[N];
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
  for (s in 1:S) {
    phi[s] = exp(log_phi[s]);
  }
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
generated quantities {
  vector[N] log_lik;
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

    log_lik[n] = dirichlet_lpdf(pseg | phi[stress_id[n]] * softmax(eta));
  }
}
'

# Ensure Stan file ends with a trailing newline
writeLines(paste0(stan_dirichlet, "\n"), "softmax_dirichlet.stan")

sm <- rstan::stan_model("softmax_dirichlet.stan")

stan_data <- list(
  N         = N,
  K         = K,
  J         = J,
  S         = S,
  start_idx = start_idx,
  len       = len,
  p_obs     = as.vector(p_obs),
  g_z       = as.vector(g_z),
  tax_of    = tax_of,
  stress_id = stress_id,
  Rdim      = Rdim,
  rich_id   = rich_id,
  use_taxon_biases = as.integer(USE_TAXON_BIASES)
)

fit <- rstan::sampling(
  sm,
  data    = stan_data,
  chains  = CHAINS,
  iter    = ITER,
  warmup  = WARMUP,
  control = ADAPT,
  seed    = 123,
  cores   = 1   # run chains sequentially -> avoids SIGPIPE / PSOCK issues
)

# ---------------------- STAN DIAGNOSTICS (publication-ready) ----------
message("[4b/6] Stan sampling diagnostics ...")

# 1) Convergence diagnostics (Rhat, n_eff) for key parameters
sum_fit <- summary(fit, pars = c("mu_kappa","sigma_kappa","sigma_delta","mu_phi","sigma_phi"),
                   probs = c(0.05, 0.5, 0.95))$summary

diag_tbl <- tibble::tibble(
  Parameter = rownames(sum_fit),
  Mean      = sum_fit[,"mean"],
  SD        = sum_fit[,"sd"],
  Q05       = sum_fit[,"5%"],
  Q50       = sum_fit[,"50%"],
  Q95       = sum_fit[,"95%"],
  n_eff     = sum_fit[,"n_eff"],
  Rhat      = sum_fit[,"Rhat"]
)

readr::write_csv(diag_tbl, P_TAB("Table_S14_stan_keyparam_diagnostics.csv"))

# 2) Sampler health (divergences, treedepth hits)
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)

n_div <- sum(vapply(sp, function(x) sum(x[,"divergent__"]), numeric(1)))
n_td  <- sum(vapply(sp, function(x) sum(x[,"treedepth__"] >= ADAPT$max_treedepth), numeric(1)))
n_it  <- sum(vapply(sp, function(x) nrow(x), integer(1)))

sampler_diag <- tibble::tibble(
  metric = c("iterations_postwarmup", "divergent_transitions", "max_treedepth_hits"),
  value  = c(n_it, n_div, n_td)
)
readr::write_csv(sampler_diag, P_TAB("Table_SXX_stan_sampler_diagnostics.csv"))

if (n_div > 0) message("WARNING: divergent transitions detected: ", n_div)
if (n_td  > 0) message("WARNING: max treedepth hits detected: ", n_td)

# 3) Quick Rhat flag (publication sanity): report any Rhat > 1.01
bad_rhat <- diag_tbl %>% dplyr::filter(is.finite(Rhat), Rhat > 1.01)
if (nrow(bad_rhat) > 0) {
  message("WARNING: Rhat > 1.01 for: ", paste(bad_rhat$Parameter, collapse = ", "))
}


# ---------------------- POST. DRAW SUMMARIES -------------------------
log_lik <- rstan::extract(fit, pars = "log_lik")$log_lik  # draws x N

loo_res <- NULL
if (!is.null(log_lik) && length(dim(log_lik)) == 2) {
  # LOO expects matrix-like: iterations x observations
  loo_res <- tryCatch(
    loo::loo(log_lik),
    error = function(e) {
      message("LOO failed: ", conditionMessage(e))
      NULL
    }
  )
}

if (!is.null(loo_res)) print(loo_res)

draws <- rstan::extract(fit)
nd    <- if (!is.null(draws$kappa)) dim(draws$kappa)[1] else length(draws$mu_kappa)

set.seed(123)
sel    <- if (nd > NDRAWS_COMPOSITION) sample.int(nd, NDRAWS_COMPOSITION) else seq_len(nd)
set.seed(123)
sel_mi <- if (nd > NDRAWS_MI)          sample.int(nd, NDRAWS_MI)          else seq_len(nd)

idx_split <- split(seq_len(J), rep.int(seq_len(N), times = len))
g_z_vec   <- as.numeric(g_z)
g_raw_vec <- as.numeric(g_raw)

get_kappa <- function(d, s, r = 1L) {
  kd <- draws$kappa
  if (length(dim(kd)) == 2L) {
    kd[d, s]
  } else if (length(dim(kd)) == 3L) {
    kd[d, s, r]
  } else {
    stop("Unexpected 'kappa' dims")
  }
}

# Posterior compositions & AWM per sample
pred_list <- vector("list", N)
AWM_list  <- vector("list", N)

for (n in seq_len(N)) {
  seg     <- idx_split[[n]]
  s       <- stress_id[n]
  r_idx   <- if (USE_RICHNESS_KAPPA) rich_id[n] else 1L
  gsub_z  <- g_z_vec[seg]
  gsub_r  <- g_raw_vec[seg]
  taxa_ix <- tax_of[seg]
  
  P <- matrix(NA_real_, nrow = length(sel), ncol = length(seg))
  A <- numeric(length(sel))
  
  for (ii in seq_along(sel)) {
    d    <- sel[ii]
    kap  <- get_kappa(d, s, r_idx)
    logits <- kap * gsub_z + draws$delta[d, taxa_ix]
    p    <- softmax(logits)
    P[ii, ] <- p
    A[ii]   <- sum(p * gsub_r)
  }
  
  p_mean <- colMeans(P)
  p_q    <- apply(P, 2, stats::quantile,
                  probs = c(0.05, 0.5, 0.95),
                  na.rm = TRUE)
  
  pred_list[[n]] <- tibble::tibble(
    SampleID  = SAMPLES$SampleID[n],
    Com_Id    = SAMPLES$Com_Id[n],
    Stress    = SAMPLES$Stress[n],
    Diversity = SAMPLES$Diversity[n],
    Id        = taxa_levels[taxa_ix],
    p_hat     = p_mean,
    p_lwr     = p_q[1, ],
    p_med     = p_q[2, ],
    p_upr     = p_q[3, ]
  )
  
  AWM_list[[n]] <- tibble::tibble(
    SampleID  = SAMPLES$SampleID[n],
    Com_Id    = SAMPLES$Com_Id[n],
    Stress    = SAMPLES$Stress[n],
    Diversity = SAMPLES$Diversity[n],
    AWM_mean  = mean(A, na.rm = TRUE),
    AWM_lwr   = stats::quantile(A, 0.05, na.rm = TRUE),
    AWM_upr   = stats::quantile(A, 0.95, na.rm = TRUE)
  )
}

pred_comp <- dplyr::bind_rows(pred_list)
AWM_post  <- dplyr::bind_rows(AWM_list)

saveRDS(pred_comp, P_RDS("bayes_pred_comp.rds"))
saveRDS(AWM_post,  P_RDS("bayes_awms.rds"))

# Join observed compositions; build composition metrics
top_obs <- ord %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(p_obs = Abundance) %>%
  dplyr::ungroup() %>%
  dplyr::select(SampleID, Id, p_obs)

pred_comp <- pred_comp %>%
  dplyr::left_join(top_obs, by = c("SampleID","Id"))

comp_plot_df <- pred_comp %>%
  dplyr::filter(!is.na(p_obs)) %>%
  dplyr::mutate(Diversity = factor(Diversity))

# ---------------------- Interval coverage diagnostic -------------------
coverage_tbl <- comp_plot_df %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(
    n = dplyr::n(),
    coverage_90 = mean(p_obs >= p_lwr & p_obs <= p_upr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(coverage_90 = round(coverage_90, 3)) %>%
  dplyr::arrange(dplyr::desc(coverage_90))

readr::write_csv(coverage_tbl, P_TAB("Table_SXX_composition_interval_coverage90_by_stress.csv"))

comp_per_sample <- comp_plot_df %>%
  dplyr::group_by(SampleID, Stress, Diversity) %>%
  dplyr::summarise(
    RMSE = sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE)),
    JS   = suppressWarnings(philentropy::JSD(rbind(p_obs, p_hat))),
    .groups = "drop"
  )

comp_by_stress <- comp_per_sample %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(
    RMSE = mean(RMSE, na.rm = TRUE),
    JS   = mean(JS,   na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------- Fig S1 (facet composition) -------------------
lab_comp <- comp_by_stress %>%
  dplyr::mutate(
    label = paste0(
      "RMSE = ", sprintf("%.3f", RMSE),
      "\nJS = ",   sprintf("%.3f", JS)
    ),
    x = -Inf, y = Inf
  )

p_comp <- ggplot(comp_plot_df, aes(x = p_hat, y = p_obs, colour = Diversity)) +
  geom_abline(linetype = 2) +
  geom_errorbarh(
    aes(xmin = p_lwr, xmax = p_upr),
    alpha = 0.30, height = 0, linewidth = 0.3
  ) +
  geom_point(alpha = 0.85, size = 1.8) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  facet_wrap(~ Stress, nrow = 2) +
  scale_colour_manual(
    values = c("2"="#56B4E9","4"="#F0E442","8"="grey40"),
    name   = "Richness"
  ) +
  labs(
    x = "Predicted relative abundance (posterior mean ± 90% CI)",
    y = "Observed relative abundance"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    panel.border     = element_rect(color="black", fill=NA, linewidth=0.4),
    panel.spacing    = grid::unit(1.5, "lines")
  ) +
  geom_text(
    data        = lab_comp,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = -0.05,
    vjust       = 1.4,
    size        = 3
  )

ggsave(
  P_FIG("Fig_S4.tiff"),
  label_plot(p_comp, "S4"),
  width  = 175,
  height = 120,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# ========================= OD predictions ============================
message("[5/6] Building OD predictors and figures ...")

comm_od <- read.csv(P_IN("CommunityOD_All_Jan23_edit_zeros.csv")) %>%
  dplyr::mutate(
    Stress = dplyr::case_when(
      Stress == "pH_Sal"      ~ "pHSal",
      Stress == "pH_Sal_Temp" ~ "pHSalTemp",
      Stress == "pH_Temp"     ~ "pHTemp",
      Stress == "Sal_Temp"    ~ "SalTemp",
      TRUE                    ~ Stress
    ),
    Com_Id    = paste0("Com", gsub("\\s+","", Community)),
    Diversity = as.character(Diversity_Level)
  ) %>%
  dplyr::filter(Evo_Treatment != "Evo", t == 2) %>%
  dplyr::transmute(
    Com_Id = gsub("\\s+","", Com_Id),
    Stress, Diversity,
    OD = as.numeric(OD)
  ) %>%
  dplyr::mutate(
    Stress   = factor(Stress),
    Diversity= factor(Diversity)
  ) %>%
  droplevels()

lvl_stress <- levels(comm_od$Stress)
lvl_div    <- levels(comm_od$Diversity)

# Use aligned 'ord' with added column to avoid recycling
ord_od <- ord %>%
  dplyr::mutate(
    Stress   = factor(Stress, levels = lvl_stress),
    Diversity= factor(as.character(Diversity), levels = lvl_div),
    g_raw_col= g_raw
  )

# (1) Observed AWM (oracle) -> S2a
AWM_obs <- ord_od %>%
  dplyr::group_by(SampleID, Com_Id, Stress, Diversity) %>%
  dplyr::summarise(
    AWM_obs = sum(Abundance * g_raw_col, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Com_Id, Stress, Diversity) %>%
  dplyr::summarise(
    AWM_obs = mean(AWM_obs, na.rm = TRUE),
    .groups = "drop"
  )

od_obs <- comm_od %>%
  dplyr::left_join(AWM_obs, by = c("Com_Id","Stress","Diversity")) %>%
  dplyr::filter(is.finite(OD), is.finite(AWM_obs))

m_obs <- lm(OD ~ AWM_obs * Stress + Diversity, data = od_obs)
od_obs$pred <- predict(m_obs, newdata = od_obs)

rmse_obs <- sqrt(mean((od_obs$OD - od_obs$pred)^2))
r2_obs   <- 1 - sum((od_obs$OD - od_obs$pred)^2) /
  sum((od_obs$OD - mean(od_obs$OD))^2)

theme_nat <- theme_classic(base_family = "Helvetica") +
  theme(
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.title   = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

shape_map <- c("2" = 16, "4" = 17, "8" = 15)  # ● ▲ ■

# --- Fig S2a (Observed AWM × g) ---
p_obs <- ggplot(od_obs, aes(x = pred, y = OD, colour = Stress, shape = Diversity)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.9) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("RMSE = %.2f\nR² = %.2f", rmse_obs, r2_obs),
    hjust = -0.05, vjust = 1.4, size = 4
  ) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    title = "Observed AWM × g",
    x     = expression(Predicted~OD[600]),
    y     = expression(OD[600])
  ) +
  theme_nat

ggsave(
  P_FIG("Fig_S5a.tiff"),
  label_plot(p_obs, "S5a"),
  width  = 170,
  height = 130,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# Helper: compute AWM by draw
compute_one_AWM_draw <- function(d) {
  AWM_draw <- vector("list", N)
  for (n in seq_len(N)) {
    seg <- idx_split[[n]]
    s   <- stress_id[n]
    r   <- if (USE_RICHNESS_KAPPA) rich_id[n] else 1L
    kap <- get_kappa(d, s, r)
    logits <- kap * g_z_vec[seg] + draws$delta[d, tax_of[seg]]
    p <- softmax(logits)
    AWM_draw[[n]] <- tibble::tibble(
      SampleID = SAMPLES$SampleID[n],
      AWM      = sum(p * g_raw_vec[seg])
    )
  }
  dplyr::bind_rows(AWM_draw) %>%
    dplyr::left_join(SAMPLES, by = "SampleID") %>%
    dplyr::group_by(Com_Id, Stress, Diversity) %>%
    dplyr::summarise(AWM = mean(AWM), .groups = "drop") %>%
    dplyr::mutate(
      Stress   = factor(Stress,   levels = lvl_stress),
      Diversity= factor(as.character(Diversity), levels = lvl_div)
    )
}

# (2) Bayes AWM via MI -> Fig 3b
awm_point <- {
  tmp <- lapply(sel_mi, compute_one_AWM_draw)
  dplyr::bind_rows(tmp, .id = "draw") %>%
    dplyr::group_by(Com_Id, Stress, Diversity) %>%
    dplyr::summarise(AWM_mean = mean(AWM), .groups = "drop")
}

base_df <- comm_od %>%
  dplyr::left_join(awm_point, by = c("Com_Id","Stress","Diversity")) %>%
  dplyr::filter(is.finite(OD), is.finite(AWM_mean))

choose_od_formula <- function(df, awm_col = "AWM_mean", k_outer = 5, k_inner = 3) {
  set.seed(123)
  formulas <- list(
    f1 = as.formula(paste("OD ~", awm_col, "* Stress + Diversity")),
    f2 = as.formula(paste("OD ~", awm_col, "+ Stress + Diversity")),
    f3 = as.formula(paste("OD ~", awm_col)),
    f4 = as.formula("OD ~ Stress + Diversity")
  )
  n    <- nrow(df)
  idx  <- sample.int(n)
  folds <- split(idx, cut(seq_along(idx), breaks = k_outer, labels = FALSE))
  
  mean_rmse <- rep(0, length(formulas))
  
  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- setdiff(idx, test_idx)
    dtrain    <- df[train_idx, , drop = FALSE]
    
    inner_idx   <- sample(seq_len(nrow(dtrain)))
    inner_folds <- split(inner_idx, cut(seq_along(inner_idx),
                                        breaks = k_inner, labels = FALSE))
    rmse_mat <- matrix(NA_real_, nrow = length(formulas), ncol = length(inner_folds))
    
    for (j in seq_along(formulas)) {
      for (ii in seq_along(inner_folds)) {
        v_idx  <- inner_folds[[ii]]
        tr_idx <- setdiff(inner_idx, v_idx)
        fit <- tryCatch(
          lm(formulas[[j]], data = dtrain[tr_idx, , drop = FALSE]),
          error = function(e) NULL
        )
        if (is.null(fit)) next
        pred <- predict(fit, newdata = dtrain[v_idx, , drop = FALSE])
        rmse_mat[j, ii] <- sqrt(mean((dtrain$OD[v_idx] - pred)^2, na.rm = TRUE))
      }
    }
    
    inner_rmse <- rowMeans(rmse_mat, na.rm = TRUE)
    best_j     <- which.min(inner_rmse)
    
    fit_best <- tryCatch(
      lm(formulas[[best_j]], data = dtrain),
      error = function(e) NULL
    )
    if (is.null(fit_best)) next
    
    pred_te <- predict(fit_best, newdata = df[test_idx, , drop = FALSE])
    mean_rmse[best_j] <- mean_rmse[best_j] +
      sqrt(mean((df$OD[test_idx] - pred_te)^2, na.rm = TRUE))
  }
  
  formulas[[which.min(mean_rmse)]]
}

best_formula <- if (USE_NESTED_CV_FOR_OD) {
  choose_od_formula(base_df, awm_col = "AWM_mean", k_outer = K_OUTER, k_inner = K_INNER)
} else {
  as.formula("OD ~ AWM_mean * Stress + Diversity")
}

mi_preds <- vector("list", length(sel_mi))
keep     <- 0L

for (ii in seq_along(sel_mi)) {
  d     <- sel_mi[ii]
  awm_d <- compute_one_AWM_draw(d)
  
  dat_d <- comm_od %>%
    dplyr::left_join(awm_d, by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::filter(is.finite(OD), is.finite(AWM)) %>%
    dplyr::rename(AWM_mean = AWM)
  
  if (nrow(dat_d) < 5) next
  
  m_d <- tryCatch(lm(best_formula, data = dat_d), error = function(e) NULL)
  if (is.null(m_d)) next
  
  keep <- keep + 1L
  mi_preds[[keep]] <- tibble::tibble(
    Com_Id   = dat_d$Com_Id,
    Stress   = dat_d$Stress,
    Diversity= dat_d$Diversity,
    OD       = dat_d$OD,
    pred     = as.numeric(predict(m_d, newdata = dat_d))
  )
}

stopifnot(keep > 0)

mi_summary <- mi_preds[seq_len(keep)] %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(Com_Id, Stress, Diversity, OD) %>%
  dplyr::summarise(
    pred_mean = mean(pred),
    pred_sd   = sd(pred),
    .groups   = "drop"
  )

rmse_bayes <- sqrt(mean((mi_summary$OD - mi_summary$pred_mean)^2))
r2_bayes   <- 1 - sum((mi_summary$OD - mi_summary$pred_mean)^2) /
  sum((mi_summary$OD - mean(mi_summary$OD))^2)

p_bayes <- ggplot(mi_summary, aes(pred_mean, OD, colour = Stress, shape = Diversity)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.9, size = 2) +
  geom_errorbarh(
    aes(xmin = pred_mean - 1.96 * pred_sd,
        xmax = pred_mean + 1.96 * pred_sd),
    height = 0
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("RMSE = %.2f\nR² = %.2f", rmse_bayes, r2_bayes),
    hjust  = -0.05, vjust = 1.4, size = 4
  ) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    title = "Bayesian AWM × g (posterior MI)",
    x     = expression(Predicted~OD[600]),
    y     = expression(OD[600])
  ) +
  theme_nat

ggsave(
  P_FIG("Fig_3b.tiff"),
  label_plot(p_bayes, "3b"),
  width  = 170,
  height = 130,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# (3) Equal-abundance baseline -> S2b
df_eq <- {
  AWM_equal <- {
    out <- vector("list", N)
    for (n in seq_len(N)) {
      seg    <- idx_split[[n]]
      AWM_eq <- mean(g_raw_vec[seg], na.rm = TRUE)
      out[[n]] <- tibble::tibble(
        SampleID = SAMPLES$SampleID[n],
        AWM_eq   = AWM_eq
      )
    }
    dplyr::bind_rows(out) %>%
      dplyr::left_join(SAMPLES, by = "SampleID") %>%
      dplyr::group_by(Com_Id, Stress, Diversity) %>%
      dplyr::summarise(
        AWM_eq = mean(AWM_eq, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        Stress   = factor(Stress, levels = lvl_stress),
        Diversity= factor(as.character(Diversity), levels = lvl_div)
      )
  }
  comm_od %>%
    dplyr::left_join(AWM_equal, by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::filter(is.finite(OD), is.finite(AWM_eq)) %>%
    dplyr::rename(AWM_mean = AWM_eq)
}

m_eq <- lm(OD ~ AWM_mean * Stress + Diversity, data = df_eq)
df2  <- df_eq %>%
  dplyr::mutate(pred = predict(m_eq, newdata = df_eq))

rmse2 <- sqrt(mean((df2$OD - df2$pred)^2))
r2_2  <- 1 - sum((df2$OD - df2$pred)^2) /
  sum((df2$OD - mean(df2$OD))^2)

p_eq <- ggplot(df2, aes(pred, OD, colour = Stress, shape = Diversity)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = .9) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("RMSE = %.2f\nR² = %.2f", rmse2, r2_2),
    hjust  = -0.05, vjust = 1.4, size = 4
  ) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    title = "Equal-abundance baseline (×g)",
    x     = expression(Predicted~OD[600]),
    y     = expression(OD[600])
  ) +
  theme_nat

ggsave(
  P_FIG("Fig_S5b.tiff"),
  label_plot(p_eq, "S5b"),
  width  = 170,
  height = 130,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# Pool of growth rates by stress (for trait-shuffled null)
g_pool_by_stress <- ord %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(g_pool = list(g), .groups = "drop")

g_pool_map <- setNames(g_pool_by_stress$g_pool,
                       as.character(g_pool_by_stress$Stress))

# (4) Trait-shuffled baseline -> S2c
df_rand <- {
  RANDOM_DRAWS <- 200  # enough to stabilise, but not huge
  
  AWM_random <- {
    out <- vector("list", N)
    set.seed(123)
    
    for (n in seq_len(N)) {
      seg <- idx_split[[n]]
      L   <- length(seg)
      if (L < 1) next
      
      # Stress for this community
      s_n    <- SAMPLES$Stress[n]
      g_pool <- g_pool_map[[as.character(s_n)]]
      
      # If for some reason pool is empty, skip
      if (length(g_pool) < 1) next
      
      # Randomly assign growth rates from the stress-specific pool
      A <- replicate(RANDOM_DRAWS, {
        g_samp <- sample(g_pool, size = L, replace = TRUE)
        mean(g_samp)
      })
      
      out[[n]] <- tibble::tibble(
        SampleID = SAMPLES$SampleID[n],
        AWM_rand = mean(A, na.rm = TRUE)
      )
    }
    
    dplyr::bind_rows(out) %>%
      dplyr::left_join(SAMPLES, by = "SampleID") %>%
      dplyr::group_by(Com_Id, Stress, Diversity) %>%
      dplyr::summarise(
        AWM_rand = mean(AWM_rand, na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        Stress   = factor(Stress, levels = lvl_stress),
        Diversity= factor(as.character(Diversity), levels = lvl_div)
      )
  }
  
  comm_od %>%
    dplyr::left_join(AWM_random, by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::filter(is.finite(OD), is.finite(AWM_rand)) %>%
    dplyr::rename(AWM_mean = AWM_rand)
}

m_rand <- lm(OD ~ AWM_mean * Stress + Diversity, data = df_rand)
df3    <- df_rand %>%
  dplyr::mutate(pred = predict(m_rand, newdata = df_rand))

rmse3 <- sqrt(mean((df3$OD - df3$pred)^2))
r2_3  <- 1 - sum((df3$OD - df3$pred)^2) /
  sum((df3$OD - mean(df3$OD))^2)

p_rand <- ggplot(df3, aes(pred, OD, colour = Stress, shape = Diversity)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = .9) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("RMSE = %.2f\nR² = %.2f", rmse3, r2_3),
    hjust  = -0.05, vjust = 1.4, size = 4
  ) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    title = "Trait-shuffled baseline (×g)",
    x     = expression(Predicted~OD[600]),
    y     = expression(OD[600])
  ) +
  theme_nat

ggsave(
  P_FIG("Fig_S5c.tiff"),
  label_plot(p_rand, "S5c"),
  width  = 170,
  height = 130,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# ---------------- Single-panel composition (Fig 3a) ------------------
w    <- comp_plot_df$p_obs
y    <- comp_plot_df$p_obs
yhat <- comp_plot_df$p_hat

wSSE   <- sum(w * (y - yhat)^2, na.rm = TRUE)
ybar_w <- sum(w * y, na.rm = TRUE) / sum(w, na.rm = TRUE)
wSST   <- sum(w * (y - ybar_w)^2, na.rm = TRUE)

wRMSE <- sqrt(wSSE / sum(w, na.rm = TRUE))
wR2   <- 1 - (wSSE / wSST)

FIG_comp_single <- ggplot(
  comp_plot_df,
  aes(x = p_hat, y = p_obs, colour = Stress, shape = Diversity)
) +
  geom_abline(linetype = 2) +
  geom_errorbarh(
    aes(xmin = p_lwr, xmax = p_upr),
    height = 0, alpha = 0.35, linewidth = 0.3
  ) +
  geom_point(alpha = 0.9, size = 1.9) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("wRMSE = %.2f\nwR² = %.2f", wRMSE, wR2),
    hjust  = -0.05, vjust = 1.4, size = 4
  ) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    x = "Predicted relative abundance (posterior mean ± 90% CI)",
    y = "Observed relative abundance"
  ) +
  theme_nat +
  theme(legend.position = "right")

ggsave(
  P_FIG("Fig_3a.tiff"),
  label_plot(FIG_comp_single, "3a"),
  width  = 170,
  height = 130,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# ------------------------ Tables S9–S13 -------------------------------
message("[6/6] Writing tables S9–S13 ...")

table_s9 <- comp_plot_df %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(
    n    = dplyr::n(),
    RMSE = sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE)),
    R2   = {
      mu <- mean(p_obs, na.rm = TRUE)
      1 - sum((p_obs - p_hat)^2, na.rm = TRUE) /
        sum((p_obs - mu)^2,     na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    RMSE = round(RMSE, 3),
    R2   = round(R2,   3)
  ) %>%
  dplyr::arrange(dplyr::desc(R2))

readr::write_csv(table_s9, P_TAB("Table_S9_composition_fit_metrics_by_stress.csv"))


if (!is.null(log_lik)) {
  sample_stress <- SAMPLES$Stress
  split_idx     <- split(seq_len(N), sample_stress)
  
  loo_by_stress <- tibble::tibble(
    Stress   = names(split_idx),
    elpd_loo = purrr::map_dbl(
      split_idx,
      function(idx) {
        ll_sub <- log_lik[, idx, drop = FALSE]
        suppressMessages(loo::loo(ll_sub))$estimates["elpd_loo","Estimate"]
      }
    ),
    n        = purrr::map_int(split_idx, length)
  ) %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~round(., 3))
    ) %>%
    dplyr::arrange(dplyr::desc(elpd_loo))
  
  readr::write_csv(loo_by_stress, P_TAB("Table_S10_loo_elpd_by_stress.csv"))
}

overall_tbl <- tibble::tibble(
  Model = c(
    "Bayes AWM (MI)",
    "Observed AWM (oracle ceiling)",
    "Equal-abundance × g",
    "Random × g"
  ),
  RMSE = c(
    sqrt(mean((mi_summary$OD - mi_summary$pred_mean)^2)),
    sqrt(mean((od_obs$OD      - od_obs$pred)^2)),
    sqrt(mean((df2$OD         - df2$pred)^2)),
    sqrt(mean((df3$OD         - df3$pred)^2))
  )
) %>%
  dplyr::mutate(
    R2 = c(
      {
        mu <- mean(mi_summary$OD)
        1 - sum((mi_summary$OD - mi_summary$pred_mean)^2) /
          sum((mi_summary$OD - mu)^2)
      },
      {
        mu <- mean(od_obs$OD)
        1 - sum((od_obs$OD - od_obs$pred)^2) /
          sum((od_obs$OD - mu)^2)
      },
      {
        mu <- mean(df2$OD)
        1 - sum((df2$OD - df2$pred)^2) /
          sum((df2$OD - mu)^2)
      },
      {
        mu <- mean(df3$OD)
        1 - sum((df3$OD - df3$pred)^2) /
          sum((df3$OD - mu)^2)
      }
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(where(is.numeric), ~round(., 3))
  )

readr::write_csv(overall_tbl, P_TAB("Table_S11_model_overall_metrics.csv"))

readr::write_csv(
  comp_by_stress %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~round(., 3))
    ) %>%
    dplyr::arrange(Stress),
  P_TAB("Table_S12_js_by_stress.csv")
)

if (!is.null(loo_res)) {
  est <- tibble::as_tibble(loo_res$estimates, rownames = "Stat") %>%
    dplyr::select(Stat, Estimate, SE) %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~round(., 3))
    )
  
  readr::write_csv(est, P_TAB("Table_S13_loo_summary_stats.csv"))
}

message("Saved: Fig_3a.tiff, Fig_3b.tiff, Fig_S4.tiff, Fig_S5a.tiff, Fig_S5b.tiff, Fig_S5c.tiff")
message("Tables: S9–S13 written to '", normalizePath(DIR_TABLES), "'")


# ------------------------------ PART B -------------------------------
# ======================================================================
# Blocked cross-validation validation for Fig 3 (out-of-sample evaluation)
#
# What is validated:
# - Composition prediction accuracy on held-out TEST folds
# - OD prediction accuracy on held-out TEST folds using TRAIN-only OD regression
#
# Blocking options:
# - BLOCK_TYPE = "Com_Id": K-fold blocked CV where whole communities (Com_Id)
#   are held out by fold (recommended default).
# - BLOCK_TYPE = "Stress": leave-one-stress-out CV (each Stress is a fold).
#
# CV figures (TEST folds only):
# - Fig S6: Composition fits by Stress using held-out test-fold predictions
# - Fig S7: OD predictions using held-out test-fold predictions
#
# Tables (Part B; numbering aligned to your intent):
# - Table S15: Blocked CV metrics by fold (composition RMSE/JS + OD RMSE/R²)
# - Table S16: Blocked CV composition metrics by Stress (averaged across test samples)
# - Table S17: OD distribution summary by fold (diagnostic; min/q25/median/q75/max/sd)
# - Table S18: Overall blocked CV summary (means ± SD across folds)
#
# Notes:
# - Output filenames include a tag:
#     tag = "blockedcv_<block_type>_bias<0/1>"
#   e.g. "blockedcv_com_id_bias1" or "blockedcv_stress_bias0"
# ======================================================================

message("============================================================")
message("STARTING BLOCKED CV (Part B)")
message("============================================================")

# --- Setup (safe to re-run) ---
source("scripts/utils_functions.R")
ensure_packages()

suppressPackageStartupMessages({
  library(rstan)
  library(readr)
  library(philentropy)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(grid)     
  library(ragg)     
})


# ---- Plot style defs (standalone; do not rely on Part A objects) ----
shape_map <- c("2" = 16, "4" = 17, "8" = 15)  # ● ▲ ■

theme_nat <- theme_classic(base_family = "Helvetica") +
  theme(
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.title   = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

rstan::rstan_options(auto_write = TRUE)
options(stringsAsFactors = FALSE)
set.seed(123)

# ------------------------ USER OPTIONS --------------------------------
BLOCK_TYPE <- "Com_Id"     # "Com_Id" (recommended) or "Stress"
K_FOLDS    <- 5            # only used if BLOCK_TYPE == "Com_Id"

USE_RICHNESS_KAPPA_CV <- FALSE
USE_TAXON_BIASES_CV   <- TRUE # run twice: TRUE then FALSE

# Stan sampling (reduce for quick test)
CHAINS_CV <- 4
ITER_CV   <- 2000
WARMUP_CV <- 1000
ADAPT_CV  <- list(adapt_delta = 0.99, max_treedepth = 15)

# Posterior draw subsampling
NDRAWS_COMP <- 250
NDRAWS_MI   <- 250

# OD regression formula (fit on TRAIN only, predict on TEST)
OD_FORMULA <- stats::as.formula("OD ~ AWM * Stress + Diversity")

# ======================================================================
# DESIGN (planned members)  (reuse from Part A if present; else rebuild)
# ======================================================================
if (!exists("design")) {
  design <- tibble::tribble(
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

# ======================================================================
# LOAD & PREP ABUNDANCE + GROWTH  (do NOT rely on Part A objects)
# ======================================================================
message("[CV 1/5] Loading abundance ...")

abund_long_cv <- utils::read.csv(P_IN("All_Abundance_Data.csv")) %>%
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
    Stress = dplyr::case_when(
      Stress == "pH_Sal"       ~ "pHSal",
      Stress == "pH_Sal_Temp"  ~ "pHSalTemp",
      Stress == "pH_Temp"      ~ "pHTemp",
      Stress == "Sal_Temp"     ~ "SalTemp",
      TRUE                     ~ Stress
    )
  ) %>%
  dplyr::filter(Diversity %in% c(2,4,8)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Abundance = Abundance / sum(Abundance)) %>%
  dplyr::ungroup()

abund_complete_cv <- abund_long_cv %>%
  dplyr::distinct(SampleID, Com_Id, Stress, Diversity) %>%
  dplyr::left_join(design, by = "Com_Id", relationship = "many-to-many") %>%
  dplyr::left_join(
    abund_long_cv,
    by = c("SampleID","Com_Id","Stress","Diversity","Id"),
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(Abundance = tidyr::replace_na(Abundance, 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Abundance = Abundance / sum(Abundance)) %>%
  dplyr::ungroup()

message("[CV 2/5] Loading growth ...")

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
    dplyr::select(-TempLevel, -Rep)
}

mut_pH <- function(d) {
  d %>%
    dplyr::mutate(
      Stress = dplyr::case_when(
        (SalinityLevel == 0  & TempLevel == 20) ~ "Control",
        (SalinityLevel == 0  & TempLevel == 38) ~ "Temp",
        (SalinityLevel == 20 & TempLevel == 20) ~ "Sal",
        (SalinityLevel == 20 & TempLevel == 38) ~ "SalTemp",
        TRUE ~ "DROP"
      )
    ) %>%
    dplyr::filter(Stress != "DROP") %>%
    dplyr::select(-SalinityLevel, -TempLevel, -X, -X.1) %>%
    dplyr::filter(pH_Level %in% c(5.5, 7.2)) %>%
    dplyr::mutate(Stress = dplyr::if_else(pH_Level == 5.5, paste0("pH", Stress), Stress)) %>%
    dplyr::mutate(Stress = dplyr::if_else(Stress == "pHControl", "pH", Stress)) %>%
    dplyr::select(-pH_Level, -Rep)
}

mut_sal <- function(d) {
  d %>%
    dplyr::mutate(
      Stress = dplyr::case_when(
        (pH_Level == 7.2 & TempLevel == 20) ~ "Control",
        (pH_Level == 7.2 & TempLevel == 38) ~ "Temp",
        (pH_Level == 5.5 & TempLevel == 20) ~ "pH",
        (pH_Level == 5.5 & TempLevel == 38) ~ "pHTemp",
        TRUE ~ "DROP"
      )
    ) %>%
    dplyr::filter(Stress != "DROP") %>%
    dplyr::select(-pH_Level, -TempLevel, -X) %>%
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
    dplyr::select(-SalinityLevel, -Rep)
}

all_growth_raw_cv <- dplyr::bind_rows(
  load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv"), mut_temp),
  load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv"),   mut_pH),
  load_growth_file(P_IN("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv"),         mut_sal)
)

all_growth_cv <- all_growth_raw_cv %>%
  dplyr::mutate(
    Stress = dplyr::case_when(
      Stress == "pH_Sal"       ~ "pHSal",
      Stress == "pH_Sal_Temp"  ~ "pHSalTemp",
      Stress == "pH_Temp"      ~ "pHTemp",
      Stress == "Temp_Sal"     ~ "SalTemp",
      TRUE                     ~ Stress
    )
  ) %>%
  dplyr::group_by(Id, Stress) %>%
  dplyr::summarise(g = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Id = trimf(Id)) %>%
  dplyr::group_by(Stress) %>%
  dplyr::mutate(g_z = as.numeric(scale(g))) %>%
  dplyr::ungroup()

abund_growth_cv <- abund_complete_cv %>%
  dplyr::left_join(all_growth_cv, by = c("Id","Stress"))

complete_ids_cv <- abund_growth_cv %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(all_g = all(!is.na(g)), .groups = "drop") %>%
  dplyr::filter(all_g) %>%
  dplyr::pull(SampleID)

ord_full <- abund_growth_cv %>%
  dplyr::filter(SampleID %in% complete_ids_cv) %>%
  dplyr::arrange(SampleID, Id)

SAMPLES_cv <- ord_full %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    Com_Id    = dplyr::first(Com_Id),
    Stress    = dplyr::first(Stress),
    Diversity = dplyr::first(Diversity),
    len       = dplyr::n(),
    .groups   = "drop"
  )

taxa_levels_cv <- ord_full %>%
  dplyr::distinct(Id) %>%
  dplyr::arrange(Id) %>%
  dplyr::pull(Id)

stress_lvls_cv <- sort(unique(SAMPLES_cv$Stress))
rich_lvls_cv   <- sort(unique(SAMPLES_cv$Diversity))

# ======================================================================
# STAN MODEL for CV (own file + own model object -> no collision)
# ======================================================================
message("[CV 3/5] Compiling Stan model for CV ...")

stan_code_cv <- '
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
  int<lower=1> start_idx[N];
  int<lower=1> len[N];
  vector[J] p_obs;
  vector[J] g_z;
  int<lower=1, upper=K> tax_of[J];
  int<lower=1, upper=S> stress_id[N];
  int<lower=1> Rdim;
  int<lower=1, upper=Rdim> rich_id[N];
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

writeLines(paste0(stan_code_cv, "\n"), "softmax_dirichlet_blockedcv.stan")
sm_cv <- rstan::stan_model("softmax_dirichlet_blockedcv.stan")

# ======================================================================
# OD DATA (for CV)
# ======================================================================
comm_od_cv <- utils::read.csv(P_IN("CommunityOD_All_Jan23_edit_zeros.csv")) %>%
  dplyr::mutate(
    Stress = dplyr::case_when(
      Stress == "pH_Sal"      ~ "pHSal",
      Stress == "pH_Sal_Temp" ~ "pHSalTemp",
      Stress == "pH_Temp"     ~ "pHTemp",
      Stress == "Sal_Temp"    ~ "SalTemp",
      TRUE                    ~ Stress
    ),
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

# ======================================================================
# HELPERS (same as your CV script)
# ======================================================================
softmax_vec <- function(x) {
  ex <- exp(x - max(x))
  ex / sum(ex)
}

build_ragged <- function(ord_df, sample_ids_keep,
                         taxa_levels, stress_lvls, rich_lvls,
                         USE_RICHNESS_KAPPA, use_taxon_biases_int) {
  
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

# ======================================================================
# BLOCKED FOLDS
# ======================================================================
message("[CV 4/5] Building blocked folds ...")

if (BLOCK_TYPE == "Stress") {
  SAMPLES_cv$fold <- match(SAMPLES_cv$Stress, stress_lvls_cv)
  K_eff <- length(stress_lvls_cv)
} else {
  set.seed(123)
  blocks <- sample(unique(SAMPLES_cv$Com_Id))
  fold_map <- setNames(rep_len(seq_len(K_FOLDS), length(blocks)), blocks)
  SAMPLES_cv$fold <- as.integer(fold_map[SAMPLES_cv$Com_Id])
  K_eff <- K_FOLDS
}

# ======================================================================
# BLOCKED CV LOOP
# ======================================================================
message("[CV 5/5] Running blocked CV: ", BLOCK_TYPE, " | taxon_bias=", USE_TAXON_BIASES_CV)

use_taxon_biases_int <- as.integer(USE_TAXON_BIASES_CV)

# ---- Collect test-fold predictions for CV figures (S3, S4) ----
cv_comp_pred_all <- vector("list", K_eff)  # test-fold composition predictions
cv_od_pred_all   <- vector("list", K_eff)  # test-fold OD predictions

fold_out <- vector("list", K_eff)

for (k in seq_len(K_eff)) {
  
  test_ids  <- SAMPLES_cv$SampleID[SAMPLES_cv$fold == k]
  train_ids <- SAMPLES_cv$SampleID[SAMPLES_cv$fold != k]
  
  message("Fold ", k, "/", K_eff, " | train=", length(train_ids), " test=", length(test_ids))
  
  train_obj <- build_ragged(ord_full, train_ids, taxa_levels_cv, stress_lvls_cv, rich_lvls_cv,
                            USE_RICHNESS_KAPPA_CV, use_taxon_biases_int)
  test_obj  <- build_ragged(ord_full, test_ids,  taxa_levels_cv, stress_lvls_cv, rich_lvls_cv,
                            USE_RICHNESS_KAPPA_CV, use_taxon_biases_int)
  
  fit_cv <- rstan::sampling(
    sm_cv,
    data    = train_obj$stan_data,
    chains  = CHAINS_CV,
    iter    = ITER_CV,
    warmup  = WARMUP_CV,
    control = ADAPT_CV,
    seed    = 500 + k,
    cores   = 1
  )
  
  draws_cv <- rstan::extract(fit_cv)
  nd <- dim(draws_cv$kappa)[1]
  
  set.seed(1000 + k)
  sel_comp <- if (nd > NDRAWS_COMP) sample.int(nd, NDRAWS_COMP) else seq_len(nd)
  set.seed(2000 + k)
  sel_mi   <- if (nd > NDRAWS_MI)   sample.int(nd, NDRAWS_MI)   else seq_len(nd)
  
  pred_test <- predict_comp_and_awm(test_obj, draws_cv, sel_comp, USE_RICHNESS_KAPPA_CV, taxa_levels_cv)
  
  # Save test-fold composition predictions for plotting (S3)
  cv_comp_pred_all[[k]] <- pred_test$pred_comp %>%
    dplyr::mutate(fold = k)
  
  comp_per_sample <- pred_test$pred_comp %>%
    dplyr::group_by(SampleID, Com_Id, Stress, Diversity) %>%
    dplyr::summarise(
      RMSE = sqrt(mean((p_obs - p_hat)^2)),
      JS   = suppressWarnings(suppressMessages(philentropy::JSD(rbind(p_obs, p_hat)))),
      .groups = "drop"
    )
  
  pred_train_mi <- predict_comp_and_awm(train_obj, draws_cv, sel_mi, USE_RICHNESS_KAPPA_CV, taxa_levels_cv)
  
  awm_train <- pred_train_mi$awm_post %>%
    dplyr::group_by(Com_Id, Stress, Diversity) %>%
    dplyr::summarise(AWM = mean(AWM_mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Stress = as.character(Stress), Diversity = as.character(Diversity))
  
  awm_test <- pred_test$awm_post %>%
    dplyr::group_by(Com_Id, Stress, Diversity) %>%
    dplyr::summarise(AWM = mean(AWM_mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Stress = as.character(Stress), Diversity = as.character(Diversity))
  
  od_train <- comm_od_cv %>%
    dplyr::left_join(awm_train, by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::filter(is.finite(OD), is.finite(AWM))
  
  od_test <- comm_od_cv %>%
    dplyr::left_join(awm_test, by = c("Com_Id","Stress","Diversity")) %>%
    dplyr::filter(is.finite(OD), is.finite(AWM))
  
  # --- IMPORTANT: lock factor levels to avoid "new levels" errors in predict() ---
  od_train <- od_train %>%
    dplyr::mutate(
      Stress    = factor(Stress,    levels = stress_lvls_cv),
      Diversity = factor(Diversity, levels = as.character(rich_lvls_cv))
    )
  
  od_test <- od_test %>%
    dplyr::mutate(
      Stress    = factor(Stress,    levels = stress_lvls_cv),
      Diversity = factor(Diversity, levels = as.character(rich_lvls_cv))
    )
  
  od_rmse <- NA_real_
  od_r2   <- NA_real_
  od_pred_df <- NULL
  
  if (nrow(od_train) >= 10 && nrow(od_test) >= 3) {
    m_od <- stats::lm(OD_FORMULA, data = od_train)
    pred_od <- stats::predict(m_od, newdata = od_test)
    
    od_rmse <- sqrt(mean((od_test$OD - pred_od)^2))
    od_r2   <- 1 - sum((od_test$OD - pred_od)^2) / sum((od_test$OD - mean(od_test$OD))^2)
    
    od_pred_df <- od_test %>%
      dplyr::mutate(
        pred_OD = as.numeric(pred_od),
        fold    = k
      )
  } else {
    od_pred_df <- od_test %>%
      dplyr::mutate(
        pred_OD = NA_real_,
        fold    = k
      )
  }
  
  # Save test-fold OD predictions for plotting (S4)
  cv_od_pred_all[[k]] <- od_pred_df
  
  
  fold_out[[k]] <- list(
    fold = k,
    comp = comp_per_sample,
    od_rmse = od_rmse,
    od_r2 = od_r2,
    n_test_samples = length(test_ids),
    n_test_od = nrow(od_test)
  )
}

# ======================================================================
# SUMMARIES + TABLES
# ======================================================================
by_fold_comp <- dplyr::bind_rows(lapply(fold_out, function(x) x$comp %>% dplyr::mutate(fold = x$fold)))

summary_by_fold <- by_fold_comp %>%
  dplyr::group_by(fold) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    comp_RMSE = mean(RMSE, na.rm = TRUE),
    comp_JS   = mean(JS,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    od_RMSE   = sapply(fold_out, `[[`, "od_rmse"),
    od_R2     = sapply(fold_out, `[[`, "od_r2"),
    n_test_od = sapply(fold_out, `[[`, "n_test_od")
  )

overall_summary <- summary_by_fold %>%
  dplyr::summarise(
    BLOCK_TYPE = BLOCK_TYPE,
    USE_TAXON_BIASES = USE_TAXON_BIASES_CV,
    K = dplyr::n(),
    comp_RMSE_mean = mean(comp_RMSE, na.rm = TRUE),
    comp_RMSE_sd   = stats::sd(comp_RMSE, na.rm = TRUE),
    comp_JS_mean   = mean(comp_JS,   na.rm = TRUE),
    comp_JS_sd     = stats::sd(comp_JS,   na.rm = TRUE),
    od_RMSE_mean   = mean(od_RMSE,   na.rm = TRUE),
    od_RMSE_sd     = stats::sd(od_RMSE,   na.rm = TRUE),
    od_R2_mean     = mean(od_R2,     na.rm = TRUE),
    od_R2_sd       = stats::sd(od_R2,     na.rm = TRUE),
    .groups = "drop"
  )

print(overall_summary)

by_stress <- by_fold_comp %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(
    n = dplyr::n(),
    RMSE = mean(RMSE, na.rm = TRUE),
    JS   = mean(JS,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n))

tag <- paste0("blockedcv_", tolower(BLOCK_TYPE), "_bias", as.integer(USE_TAXON_BIASES_CV))

readr::write_csv(overall_summary, P_TAB(paste0("Table_S18_", tag, "_summary.csv")))
readr::write_csv(summary_by_fold,  P_TAB(paste0("Table_S15_", tag, "_by_fold.csv")))
readr::write_csv(by_stress,        P_TAB(paste0("Table_S16_", tag, "_by_stress.csv")))

message("DONE. Wrote blocked CV tables with tag: ", tag)
print(overall_summary)

# ======================================================================
# FIG S3 + FIG S4: OUT-OF-SAMPLE (TEST-FOLD) PLOTS
# ======================================================================

cv_comp_pred_df <- dplyr::bind_rows(cv_comp_pred_all) %>%
  dplyr::filter(is.finite(p_obs), is.finite(p_hat))

cv_od_pred_df <- dplyr::bind_rows(cv_od_pred_all) %>%
  dplyr::filter(is.finite(OD), is.finite(pred_OD))

# ---- Fig S3: CV test composition fits by stress ----
# Correct CV composition metrics: compute per-sample vector metrics, then average by stress
cv_comp_per_sample <- cv_comp_pred_df %>%
  dplyr::group_by(SampleID, Com_Id, Stress, Diversity) %>%
  dplyr::summarise(
    RMSE = sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE)),
    JS   = suppressWarnings(suppressMessages(philentropy::JSD(rbind(p_obs, p_hat)))),
    .groups = "drop"
  )

cv_comp_by_stress <- cv_comp_per_sample %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(
    n    = dplyr::n(),
    RMSE = mean(RMSE, na.rm = TRUE),
    JS   = mean(JS,   na.rm = TRUE),
    .groups = "drop"
  )

lab_s3 <- cv_comp_by_stress %>%
  dplyr::mutate(
    label = paste0("RMSE = ", sprintf("%.3f", RMSE), "\nJS = ", sprintf("%.3f", JS)),
    x = -Inf, y = Inf
  )

p_s3 <- ggplot(cv_comp_pred_df, aes(x = p_hat, y = p_obs, colour = factor(Diversity))) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.85, size = 1.7) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  facet_wrap(~ Stress, nrow = 2) +
  scale_colour_manual(values = c("2"="#56B4E9","4"="#F0E442","8"="grey40"),
                      name = "Richness") +
  labs(
    x = "Predicted relative abundance (CV test; posterior mean)",
    y = "Observed relative abundance"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    panel.border     = element_rect(color="black", fill=NA, linewidth=0.4),
    panel.spacing    = grid::unit(1.5, "lines")
  ) +
  geom_text(
    data = lab_s3,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.4, size = 3
  )

ggsave(
  P_FIG("Fig_S6.tiff"),
  label_plot(p_s3, "S6"),
  width  = 175,
  height = 120,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# ---- Fig S4: CV test OD predictions ----
rmse_s4 <- sqrt(mean((cv_od_pred_df$OD - cv_od_pred_df$pred_OD)^2, na.rm = TRUE))
r2_s4   <- 1 - sum((cv_od_pred_df$OD - cv_od_pred_df$pred_OD)^2, na.rm = TRUE) /
  sum((cv_od_pred_df$OD - mean(cv_od_pred_df$OD, na.rm = TRUE))^2, na.rm = TRUE)

p_s4 <- ggplot(cv_od_pred_df, aes(x = pred_OD, y = OD, colour = Stress, shape = factor(Diversity))) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.85, size = 2) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("RMSE = %.2f\nR² = %.2f", rmse_s4, r2_s4),
    hjust = -0.05, vjust = 1.4, size = 4
  ) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    x = expression(Predicted~OD[600]~"(CV test)"),
    y = expression(Observed~OD[600])
  ) +
  theme_nat

ggsave(
  P_FIG("Fig_S7.tiff"),
  label_plot(p_s4, "S7"),
  width  = 170,
  height = 130,
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

# ======================================================================
# EXTRA: OD spread by fold (your diagnostic)
# ======================================================================
message("============================================================")
message("OD distribution by fold diagnostic")
message("============================================================")

# Recreate fold map exactly (Com_Id block)
if (BLOCK_TYPE == "Com_Id") {
  set.seed(123)
  blocks <- sample(unique(SAMPLES_cv$Com_Id))
  fold_map <- setNames(rep_len(seq_len(K_FOLDS), length(blocks)), blocks)
  SAMPLES_cv$fold <- as.integer(fold_map[SAMPLES_cv$Com_Id])
}

od_fold <- comm_od_cv %>%
  dplyr::left_join(
    SAMPLES_cv %>% dplyr::distinct(Com_Id, fold),
    by = "Com_Id"
  ) %>%
  dplyr::filter(!is.na(fold))

od_fold_summary <- od_fold %>%
  dplyr::group_by(fold) %>%
  dplyr::summarise(
    n = dplyr::n(),
    OD_min = min(OD, na.rm = TRUE),
    OD_q25 = quantile(OD, 0.25, na.rm = TRUE),
    OD_med = median(OD, na.rm = TRUE),
    OD_q75 = quantile(OD, 0.75, na.rm = TRUE),
    OD_max = max(OD, na.rm = TRUE),
    OD_sd  = sd(OD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(fold)

print(od_fold_summary, n = Inf)
readr::write_csv(od_fold_summary, P_TAB(paste0("Table_S17_", tag, "_od_by_fold_summary.csv")))

message("============================================================")
message("ALL DONE: Fig 3 outputs + blocked CV outputs written.")
message("============================================================")
