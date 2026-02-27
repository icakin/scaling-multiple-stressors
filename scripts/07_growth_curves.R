### Logistic/Growth models for Temp multi-stress experiment (clean pipeline)
### - fits: baranyi, baranyi_without_lag, gompertz, buchanan, buchanan_without_lag
### - selects best per curve via AICc
### - saves per-ID RDS + plots + parameter tables
### - SAFE against broom::tidy() chol2inv crash
### - OPTIONAL: supplementary multi-page PDF (one page per ID)
### - ADD-ON: single supplementary figure (all IDs; only 3 temps per Id x Stress)
###          BUT legend shows NUMERIC temperatures (e.g., 20°C, 30°C), not low/opt/high
### - Time shown in minutes on plots (set multiplier below if your raw t is not minutes)

rm(list = ls())

suppressPackageStartupMessages({
  library(nls.multstart)
  library(ggplot2)
  library(broom)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(MuMIn)     # AICc
  library(readr)
  library(tibble)
})

# -------------------------
# 0) Project-root paths + settings
# -------------------------
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_path <- file.path(ROOT, "data", "Cut_OD_data", "cutPoint_D17_edit2_temp.csv")
model_src <- file.path(ROOT, "scripts", "06_growth_curve_models.R")

out_rds_dir   <- file.path(ROOT, "results", "rds")
out_plot_dir  <- file.path(ROOT, "results", "figures")
out_table_dir <- file.path(ROOT, "results", "tables")

dir.create(out_rds_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_plot_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_table_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_path)) stop("data_path not found:\n  ", data_path)
if (!file.exists(model_src)) stop("model_src not found:\n  ", model_src)

# IMPORTANT: use YOUR original model definitions
# Must define: baranyi(), baranyi_without_lag(), gompertz(), buchanan(), buchanan_without_lag()
source(model_src)

# ---- TIME IN MINUTES ----
# If your raw column t is already minutes, keep this as 1.
# If your raw column t is hours, set time_multiplier <- 60.
time_multiplier <- 1
time_label <- "Time (min)"

# -------------------------
# 1) Load + clean data
# -------------------------
d_raw <- readr::read_csv(data_path, show_col_types = FALSE)

req_cols <- c("OD", "t", "TempLevel", "Stress", "Rep", "Id")
missing_cols <- setdiff(req_cols, names(d_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CSV: ", paste(missing_cols, collapse = ", "))
}

d <- d_raw %>%
  transmute(
    t = as.numeric(t) * time_multiplier,   # now in minutes (if multiplier set correctly)
    OD = as.numeric(OD),
    Rep = as.factor(Rep),
    TempLevel = as.numeric(TempLevel),
    Stress = as.factor(Stress),
    Id = as.character(Id),
    
    # enforce lower OD bound (avoid log10 issues)
    od_cor2 = pmax(OD, 0.040),
    log10_od_cor = log10(od_cor2)
  ) %>%
  filter(
    !is.na(log10_od_cor), is.finite(log10_od_cor),
    !is.na(t),            is.finite(t),
    !is.na(TempLevel),    is.finite(TempLevel),
    !is.na(Stress),
    !is.na(Rep),
    !is.na(Id),
    Id != "BLANK"
  )

# -------------------------
# 2) Model fitting helpers (UNCHANGED fitting behaviour)
# -------------------------
safe_nls <- function(formula, data, iter, start_lower, start_upper, lower, upper) {
  tryCatch(
    nls_multstart(
      formula = formula,
      data = data,
      iter = iter,
      start_lower = start_lower,
      start_upper = start_upper,
      supp_errors = "Y",
      na.action = na.omit,
      lower = lower,
      upper = upper
    ),
    error = function(e) NULL
  )
}

fit_all_models_one_curve <- function(df_curve) {
  
  if (nrow(df_curve) < 8) {
    return(list(
      fits_baranyi = NULL,
      fits_baranyi_without_lag = NULL,
      fits_gompertz = NULL,
      fits_buchanan = NULL,
      fits_buchanan_without_lag = NULL
    ))
  }
  
  list(
    fits_baranyi = safe_nls(
      log10_od_cor ~ baranyi(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_baranyi_without_lag = safe_nls(
      log10_od_cor ~ baranyi_without_lag(log10_nmax, log10_n0, mumax, t = t),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10)
    ),
    
    fits_gompertz = safe_nls(
      log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_buchanan = safe_nls(
      log10_od_cor ~ buchanan(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_buchanan_without_lag = safe_nls(
      log10_od_cor ~ buchanan_without_lag(log10_nmax, log10_n0, mumax, t = t),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10)
    )
  )
}

aicc_or_na <- function(model_obj) {
  if (is.null(model_obj)) return(NA_real_)
  tryCatch(MuMIn::AICc(model_obj), error = function(e) NA_real_)
}

pick_best_model <- function(model_list) {
  aiccs <- map_dbl(model_list, aicc_or_na)
  if (all(is.na(aiccs))) {
    return(list(best_name = NA_character_, best_model = NULL, best_aicc = NA_real_))
  }
  best_name <- names(which.min(aiccs))
  list(
    best_name = best_name,
    best_model = model_list[[best_name]],
    best_aicc = min(aiccs, na.rm = TRUE)
  )
}

make_pred_grid <- function(df, n = 120) {
  data.frame(t = seq(min(df$t, na.rm = TRUE), max(df$t, na.rm = TRUE), length.out = n))
}

augment_safe <- function(model_obj, newdata) {
  if (is.null(model_obj)) return(NULL)
  tryCatch(broom::augment(model_obj, newdata = newdata), error = function(e) NULL)
}

# ---- SAFE tidy (kept identical behaviour to your “safe” version)
tidy_safe <- function(model_obj) {
  if (is.null(model_obj)) return(NULL)
  
  td <- tryCatch(broom::tidy(model_obj), error = function(e) NULL)
  if (!is.null(td)) return(td)
  
  cf <- tryCatch(coef(model_obj), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  
  tibble(
    term = names(cf),
    estimate = as.numeric(cf),
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_
  )
}

# -------------------------
# 3) Run per-ID
# -------------------------
ids_to_run <- c("I9","I11","I15","I20","I22","I23","D11","D14","D17","I2","I8")

exclude_rules <- tibble::tribble(
  ~Id,  ~TempLevel, ~Stress, ~Rep
  # leave empty unless needed
)

apply_exclusions <- function(df, rules) {
  if (nrow(rules) == 0) return(df)
  df %>%
    mutate(Stress = as.character(Stress), Rep = as.character(Rep)) %>%
    anti_join(
      rules %>% mutate(Stress = as.character(Stress), Rep = as.character(Rep)),
      by = c("Id", "TempLevel", "Stress", "Rep")
    )
}

d_run <- d %>%
  filter(Id %in% ids_to_run) %>%
  apply_exclusions(exclude_rules)

all_summaries <- list()

# OPTIONAL: multi-page supplementary PDF (one page per ID)
make_supplement_pdf <- TRUE
supp_pdf_file <- file.path(out_plot_dir, "SUPP_Monoculture_CurveFits_bestModel_ALL_IDs.pdf")
if (make_supplement_pdf) {
  grDevices::pdf(supp_pdf_file, width = 12, height = 10, onefile = TRUE)
}

for (this_id in ids_to_run) {
  
  message("----- Running ID: ", this_id, " -----")
  
  df_id <- d_run %>% filter(Id == this_id)
  
  if (nrow(df_id) == 0) {
    message("No data for ", this_id, " (skipping).")
    next
  }
  
  # Fit per curve: Id × Rep × TempLevel × Stress
  models_tbl <- df_id %>%
    group_by(Id, Rep, TempLevel, Stress) %>%
    nest() %>%
    mutate(
      fits = map(data, fit_all_models_one_curve),
      best = map(fits, pick_best_model),
      best_model_name = map_chr(best, "best_name"),
      best_aicc = map_dbl(best, "best_aicc"),
      best_model = map(best, "best_model")
    ) %>%
    ungroup()
  
  rds_file <- file.path(out_rds_dir, paste0(this_id, "_LogisticMod.rds"))
  saveRDS(models_tbl, rds_file)
  message("Saved: ", rds_file)
  
  # AICc histogram
  model_stack <- models_tbl %>%
    select(Id, Rep, TempLevel, Stress, fits) %>%
    mutate(
      fits_baranyi = map(fits, "fits_baranyi"),
      fits_baranyi_without_lag = map(fits, "fits_baranyi_without_lag"),
      fits_buchanan = map(fits, "fits_buchanan"),
      fits_buchanan_without_lag = map(fits, "fits_buchanan_without_lag"),
      fits_gompertz = map(fits, "fits_gompertz")
    ) %>%
    select(-fits) %>%
    pivot_longer(
      cols = starts_with("fits_"),
      names_to = "model",
      values_to = "model_obj"
    ) %>%
    mutate(
      aic = map_dbl(model_obj, aicc_or_na),
      model_name = case_when(
        model == "fits_baranyi" ~ "(a) baranyi",
        model == "fits_baranyi_without_lag" ~ "(b) baranyi without lag",
        model == "fits_buchanan" ~ "(c) buchanan",
        model == "fits_buchanan_without_lag" ~ "(d) buchanan without lag",
        model == "fits_gompertz" ~ "(e) gompertz",
        TRUE ~ model
      )
    )
  
  aic_means <- model_stack %>%
    filter(!is.na(aic)) %>%
    group_by(model, model_name) %>%
    summarise(
      mean_aic = mean(aic),
      median_aic = median(aic),
      n = n(),
      .groups = "drop"
    )
  
  p_aic <- ggplot(model_stack, aes(aic)) +
    geom_histogram(fill = "lightgrey", col = "black", binwidth = 10) +
    geom_vline(data = aic_means, aes(xintercept = mean_aic), col = "red") +
    geom_vline(data = aic_means, aes(xintercept = median_aic), col = "blue") +
    facet_wrap(~ model_name, ncol = 2) +
    theme_bw(base_size = 14) +
    ylab("Count") + xlab("AICc score") +
    theme(strip.background = element_blank(),
          strip.text = element_text(hjust = 0))
  
  ggsave(
    filename = file.path(out_plot_dir, paste0(this_id, "_AICc_hist.pdf")),
    plot = p_aic, width = 10, height = 8
  )
  
  # Predictions using BEST model per curve
  preds_tbl <- models_tbl %>%
    mutate(
      pred_grid = map(data, make_pred_grid, n = 120),
      preds = map2(best_model, pred_grid, augment_safe)
    ) %>%
    select(Id, Rep, TempLevel, Stress, preds) %>%
    unnest(preds)
  
  # Curve plot per ID
  p_fit <- ggplot() +
    geom_point(data = df_id, aes(t, log10_od_cor, col = as.factor(Rep)), alpha = 0.8) +
    geom_line(data = preds_tbl, aes(t, .fitted, col = as.factor(Rep)), linewidth = 0.8) +
    facet_grid(TempLevel ~ Stress, scales = "free") +
    theme_bw(base_size = 12) +
    labs(title = paste0(this_id, " (best model per curve)"),
         x = time_label, y = "log10(OD)")
  
  ggsave(
    filename = file.path(out_plot_dir, paste0(this_id, "_CurveFits_bestModel.pdf")),
    plot = p_fit, width = 12, height = 10
  )
  
  if (make_supplement_pdf) print(p_fit)
  
  # Parameter table (SAFE)
  params_tbl <- models_tbl %>%
    mutate(params = map(best_model, tidy_safe)) %>%
    select(Id, Rep, TempLevel, Stress, best_model_name, best_aicc, params) %>%
    unnest(params)
  
  write.csv(
    params_tbl,
    file = file.path(out_table_dir, paste0(this_id, "_BestModel_Params.csv")),
    row.names = FALSE
  )
  
  best_counts <- models_tbl %>%
    count(best_model_name) %>%
    mutate(
      perc = round(100 * n / sum(n), 1),
      Id = this_id
    )
  
  all_summaries[[this_id]] <- best_counts
}

if (make_supplement_pdf) {
  grDevices::dev.off()
  message("Saved supplementary curve-fit PDF: ", supp_pdf_file)
}

if (length(all_summaries) > 0) {
  winner_summary <- bind_rows(all_summaries) %>%
    arrange(Id, desc(n))
  
  write.csv(
    winner_summary,
    file = file.path(out_table_dir, "BestModel_WinnerSummary_allIDs.csv"),
    row.names = FALSE
  )
}

message("Main loop done. Outputs in:\n  ", file.path(ROOT, "results"))

# =========================================================
# ADD-ON: single supplementary figure (ALL IDs in ONE PDF)
# - selects only 3 temps per Id x Stress: low / opt / high (by median mumax)
# - legend shows NUMERIC temps (e.g., 20°C, 30°C), not low/opt/high
# =========================================================
make_single_supp_allIDs <- TRUE

if (make_single_supp_allIDs) {
  
  message("Building single supplementary figure: 3 selected temps per ID x Stress (numeric legend) ...")
  
  rds_files <- file.path(out_rds_dir, paste0(ids_to_run, "_LogisticMod.rds"))
  rds_files <- rds_files[file.exists(rds_files)]
  
  if (length(rds_files) == 0) {
    warning("No RDS model files found in: ", out_rds_dir,
            "\nRun the main loop first so *_LogisticMod.rds files exist.")
  } else {
    
    all_models <- purrr::map_dfr(rds_files, readRDS)
    
    all_preds <- all_models %>%
      mutate(
        pred_grid = purrr::map(data, make_pred_grid, n = 120),
        preds = purrr::map2(best_model, pred_grid, augment_safe)
      ) %>%
      dplyr::select(Id, Rep, TempLevel, Stress, preds) %>%
      tidyr::unnest(preds) %>%
      dplyr::mutate(
        Rep = as.factor(Rep),
        Stress = as.factor(Stress),
        TempLevel = as.numeric(TempLevel)
      )
    
    # mumax per curve -> used to pick opt temperature per Id x Stress
    mumax_tbl <- all_models %>%
      mutate(params = purrr::map(best_model, tidy_safe)) %>%
      dplyr::select(Id, Rep, TempLevel, Stress, params) %>%
      tidyr::unnest(params) %>%
      dplyr::filter(term == "mumax") %>%
      dplyr::transmute(
        Id = as.character(Id),
        Rep = as.character(Rep),
        TempLevel = as.numeric(TempLevel),
        Stress = as.character(Stress),
        mumax = as.numeric(estimate)
      ) %>%
      dplyr::filter(is.finite(mumax))
    
    if (nrow(mumax_tbl) == 0) {
      warning("No mumax estimates found (term == 'mumax'). Cannot pick opt temps.")
    } else {
      
      # pick low/opt/high temps (opt = max median mumax), then make numeric label
      temp_pick <- mumax_tbl %>%
        dplyr::group_by(Id, Stress, TempLevel) %>%
        dplyr::summarise(med_mumax = median(mumax, na.rm = TRUE), .groups = "drop") %>%
        dplyr::group_by(Id, Stress) %>%
        dplyr::summarise(
          low_temp  = min(TempLevel, na.rm = TRUE),
          high_temp = max(TempLevel, na.rm = TRUE),
          opt_temp  = TempLevel[which.max(med_mumax)],
          .groups = "drop"
        ) %>%
        tidyr::pivot_longer(
          cols = c(low_temp, opt_temp, high_temp),
          names_to = "temp_cat",
          values_to = "TempLevel"
        ) %>%
        dplyr::group_by(Id, Stress, TempLevel) %>%
        dplyr::summarise(
          TempLabel = paste0(unique(TempLevel), "\u00B0C"),
          .groups = "drop"
        )
      
      d_plot <- d_run %>%
        dplyr::mutate(
          Id = as.character(Id),
          Stress = as.character(Stress),
          Rep = as.factor(Rep),
          TempLevel = as.numeric(TempLevel)
        ) %>%
        dplyr::inner_join(temp_pick, by = c("Id", "Stress", "TempLevel"))
      
      preds_plot <- all_preds %>%
        dplyr::mutate(
          Id = as.character(Id),
          Stress = as.character(Stress),
          TempLevel = as.numeric(TempLevel)
        ) %>%
        dplyr::inner_join(temp_pick, by = c("Id", "Stress", "TempLevel"))
      
      if (nrow(d_plot) == 0 || nrow(preds_plot) == 0) {
        warning("After filtering to selected temps, nothing left to plot.")
      } else {
        
        preds_plot$TempLabel <- factor(preds_plot$TempLabel, levels = sort(unique(preds_plot$TempLabel)))
        d_plot$TempLabel     <- factor(d_plot$TempLabel,     levels = levels(preds_plot$TempLabel))
        
        p_supp_all <- ggplot() +
          geom_point(
            data = d_plot,
            aes(x = t, y = log10_od_cor, colour = Rep),
            alpha = 0.45, size = 0.7
          ) +
          geom_line(
            data = preds_plot,
            aes(x = t, y = .fitted, colour = Rep, linetype = TempLabel),
            linewidth = 0.6
          ) +
          facet_grid(
            rows = vars(Id),
            cols = vars(Stress),
            scales = "free"
          ) +
          theme_bw(base_size = 10) +
          labs(
            title = "Monoculture growth curves: best-model fits (all IDs) - selected temperatures per stress",
            x = time_label,
            y = "log10(OD)",
            linetype = "Temperature",
            colour = "Rep"
          ) +
          theme(
            strip.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "right"
          )
        
        supp_single_file <- file.path(out_plot_dir, "Fig_S1.pdf")

        ggsave(
          filename = supp_single_file,
          plot = p_supp_all,
          width = 14,
          height = max(6, 0.9 * length(ids_to_run) + 3)
        )
        ggsave(
          filename = file.path(out_plot_dir, "Fig_S1.png"),
          plot = p_supp_all,
          width = 14,
          height = max(6, 0.9 * length(ids_to_run) + 3),
          dpi = 300
        )

        message("Saved single supplementary figure: ", supp_single_file)
      }
    }
  }
}


### Logistic/Growth models for Salinity experiment (project-root clean pipeline)
### - SAME curve-fit logic as your TEMP script:
###   baranyi, baranyi_without_lag, gompertz, buchanan, buchanan_without_lag
###   best per curve via AICc
### - Uses scripts/growth_curve_models.R (your existing file)
### - Project-root paths:
###     data/...      for inputs
###     results/...   for outputs
### - Reproducible per-curve starts (seeded by Id×Rep×Salinity×Stress)
### - SUPPLEMENTARY: Fig_S2.pdf (ALL IDs; low/opt/high salinity per Id×Stress)
###   legend shows NUMERIC salinity values (not low/opt/high)
### - Multi-page supplementary PDF is OFF

rm(list = ls())

suppressPackageStartupMessages({
  library(nls.multstart)
  library(ggplot2)
  library(broom)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(MuMIn)     # AICc
  library(readr)
  library(tibble)
})

# -------------------------
# 0) Project-root paths
# -------------------------
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

P_IN  <- function(...) file.path(ROOT, "data", ...)
P_SRC <- function(...) file.path(ROOT, "scripts", ...)
P_OUT <- function(...) file.path(ROOT, "results", ...)

DIR_RDS <- P_OUT("rds")
DIR_FIG <- P_OUT("figures")
DIR_TAB <- P_OUT("tables")

dir.create(DIR_RDS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 0b) Input + model source
# -------------------------
data_path <- P_IN("Cut_OD_data", "cutPoint_D17_edit3_Sal.csv")
model_src <- P_SRC("06_growth_curve_models.R")

if (!file.exists(data_path)) {
  stop("data_path not found:\n  ", data_path,
       "\n\nExpected:\n  data/Cut_OD_data/cutPoint_D17_edit3_Sal.csv",
       "\nCurrent getwd():\n  ", getwd())
}
if (!file.exists(model_src)) {
  stop("model_src not found:\n  ", model_src,
       "\n\nExpected:\n  scripts/growth_curve_models.R",
       "\nCurrent getwd():\n  ", getwd())
}

source(model_src)

# ---- TIME IN MINUTES ----
# If your raw column t is already minutes, keep this as 1.
# If your raw column t is hours, set time_multiplier <- 60.
time_multiplier <- 1
time_label <- "Time (min)"

# -------------------------
# 1) Load + clean data
#    CSV columns include: t, OD, SalinityLevel, pH_Level, Rep, Id
#    We map pH_Level -> Stress for compatibility (same as your salinity script).
# -------------------------
d_raw <- readr::read_csv(data_path, show_col_types = FALSE)

req_cols <- c("OD", "t", "SalinityLevel", "pH_Level", "Rep", "Id")
missing_cols <- setdiff(req_cols, names(d_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CSV: ", paste(missing_cols, collapse = ", "))
}

d <- d_raw %>%
  transmute(
    t = as.numeric(t) * time_multiplier,
    OD = as.numeric(OD),
    SalinityLevel = as.numeric(SalinityLevel),
    Stress = as.factor(pH_Level),
    Rep = as.factor(Rep),
    Id = as.character(Id),
    
    od_cor2 = pmax(OD, 0.040),
    log10_od_cor = log10(od_cor2)
  ) %>%
  filter(
    !is.na(log10_od_cor), is.finite(log10_od_cor),
    !is.na(t),            is.finite(t),
    !is.na(SalinityLevel), is.finite(SalinityLevel),
    !is.na(Stress),
    !is.na(Rep),
    !is.na(Id),
    Id != "BLANK"
  ) %>%
  arrange(Id, Stress, Rep, SalinityLevel, t)

# -------------------------
# 2) Deterministic seeding per curve (Id×Rep×Salinity×Stress)
#    Keeps fits stable across runs without changing the model logic.
# -------------------------
seed_from_group <- function(Id, Rep, SalinityLevel, Stress) {
  x <- paste(Id, Rep, SalinityLevel, Stress, sep = "|")
  as.integer(sum(utf8ToInt(x)) %% .Machine$integer.max)
}

# -------------------------
# 3) Model fitting helpers (SAME bounds / iter as your reference script)
# -------------------------
safe_nls <- function(formula, data, iter, start_lower, start_upper, lower, upper) {
  tryCatch(
    nls_multstart(
      formula = formula,
      data = data,
      iter = iter,
      start_lower = start_lower,
      start_upper = start_upper,
      supp_errors = "Y",
      na.action = na.omit,
      lower = lower,
      upper = upper
    ),
    error = function(e) NULL
  )
}

fit_all_models_one_curve <- function(df_curve, Id, Rep, SalinityLevel, Stress) {
  
  if (nrow(df_curve) < 8) {
    return(list(
      fits_baranyi = NULL,
      fits_baranyi_without_lag = NULL,
      fits_gompertz = NULL,
      fits_buchanan = NULL,
      fits_buchanan_without_lag = NULL
    ))
  }
  
  set.seed(seed_from_group(Id, Rep, SalinityLevel, Stress))
  
  list(
    fits_baranyi = safe_nls(
      log10_od_cor ~ baranyi(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_baranyi_without_lag = safe_nls(
      log10_od_cor ~ baranyi_without_lag(log10_nmax, log10_n0, mumax, t = t),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10)
    ),
    
    fits_gompertz = safe_nls(
      log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_buchanan = safe_nls(
      log10_od_cor ~ buchanan(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_buchanan_without_lag = safe_nls(
      log10_od_cor ~ buchanan_without_lag(log10_nmax, log10_n0, mumax, t = t),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10)
    )
  )
}

aicc_or_na <- function(model_obj) {
  if (is.null(model_obj)) return(NA_real_)
  tryCatch(MuMIn::AICc(model_obj), error = function(e) NA_real_)
}

pick_best_model <- function(model_list) {
  aiccs <- purrr::map_dbl(model_list, aicc_or_na)
  if (all(is.na(aiccs))) {
    return(list(best_name = NA_character_, best_model = NULL, best_aicc = NA_real_))
  }
  best_name <- names(which.min(aiccs))
  list(
    best_name = best_name,
    best_model = model_list[[best_name]],
    best_aicc = min(aiccs, na.rm = TRUE)
  )
}

make_pred_grid <- function(df, n = 120) {
  data.frame(t = seq(min(df$t, na.rm = TRUE), max(df$t, na.rm = TRUE), length.out = n))
}

augment_safe <- function(model_obj, newdata) {
  if (is.null(model_obj)) return(NULL)
  tryCatch(broom::augment(model_obj, newdata = newdata), error = function(e) NULL)
}

# SAFE tidy (avoids chol2inv crash)
tidy_safe <- function(model_obj) {
  if (is.null(model_obj)) return(NULL)
  
  td <- tryCatch(broom::tidy(model_obj), error = function(e) NULL)
  if (!is.null(td)) return(td)
  
  cf <- tryCatch(coef(model_obj), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  
  tibble(
    term = names(cf),
    estimate = as.numeric(cf),
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_
  )
}

# -------------------------
# 4) Run per-ID (multi-page supplementary PDF OFF)
# -------------------------
ids_to_run <- unique(d$Id)

exclude_rules <- tibble::tribble(
  ~Id, ~SalinityLevel, ~Stress, ~Rep
  # leave empty unless needed
)

apply_exclusions <- function(df, rules) {
  if (nrow(rules) == 0) return(df)
  df %>%
    mutate(Stress = as.character(Stress), Rep = as.character(Rep)) %>%
    anti_join(
      rules %>% mutate(Stress = as.character(Stress), Rep = as.character(Rep)),
      by = c("Id", "SalinityLevel", "Stress", "Rep")
    )
}

d_run <- d %>% apply_exclusions(exclude_rules)

all_summaries <- list()

for (this_id in ids_to_run) {
  
  message("----- Running ID: ", this_id, " -----")
  
  df_id <- d_run %>% filter(Id == this_id)
  
  if (nrow(df_id) == 0) {
    message("No data for ", this_id, " (skipping).")
    next
  }
  
  # Fit per curve: Id × Rep × SalinityLevel × Stress
  models_tbl <- df_id %>%
    group_by(Id, Rep, SalinityLevel, Stress) %>%
    nest() %>%
    mutate(
      fits = pmap(list(data, Id, Rep, SalinityLevel, Stress),
                  ~ fit_all_models_one_curve(..1, ..2, ..3, ..4, ..5)),
      best = map(fits, pick_best_model),
      best_model_name = map_chr(best, "best_name"),
      best_aicc = map_dbl(best, "best_aicc"),
      best_model = map(best, "best_model")
    ) %>%
    ungroup()
  
  # Save per-ID RDS
  rds_file <- file.path(DIR_RDS, paste0(this_id, "_LogisticMod_Salinity.rds"))
  saveRDS(models_tbl, rds_file)
  message("Saved: ", rds_file)
  
  # AICc histogram (per model)
  model_stack <- models_tbl %>%
    select(Id, Rep, SalinityLevel, Stress, fits) %>%
    mutate(
      fits_baranyi = map(fits, "fits_baranyi"),
      fits_baranyi_without_lag = map(fits, "fits_baranyi_without_lag"),
      fits_buchanan = map(fits, "fits_buchanan"),
      fits_buchanan_without_lag = map(fits, "fits_buchanan_without_lag"),
      fits_gompertz = map(fits, "fits_gompertz")
    ) %>%
    select(-fits) %>%
    pivot_longer(
      cols = starts_with("fits_"),
      names_to = "model",
      values_to = "model_obj"
    ) %>%
    mutate(
      aic = map_dbl(model_obj, aicc_or_na),
      model_name = case_when(
        model == "fits_baranyi" ~ "(a) baranyi",
        model == "fits_baranyi_without_lag" ~ "(b) baranyi without lag",
        model == "fits_buchanan" ~ "(c) buchanan",
        model == "fits_buchanan_without_lag" ~ "(d) buchanan without lag",
        model == "fits_gompertz" ~ "(e) gompertz",
        TRUE ~ model
      )
    )
  
  aic_means <- model_stack %>%
    filter(!is.na(aic)) %>%
    group_by(model, model_name) %>%
    summarise(
      mean_aic = mean(aic),
      median_aic = median(aic),
      n = n(),
      .groups = "drop"
    )
  
  p_aic <- ggplot(model_stack, aes(aic)) +
    geom_histogram(fill = "lightgrey", col = "black", binwidth = 10) +
    geom_vline(data = aic_means, aes(xintercept = mean_aic), col = "red") +
    geom_vline(data = aic_means, aes(xintercept = median_aic), col = "blue") +
    facet_wrap(~ model_name, ncol = 2) +
    theme_bw(base_size = 14) +
    ylab("Count") + xlab("AICc score") +
    theme(strip.background = element_blank(),
          strip.text = element_text(hjust = 0))
  
  ggsave(
    filename = file.path(DIR_FIG, paste0(this_id, "_AICc_hist_Salinity.pdf")),
    plot = p_aic, width = 10, height = 8
  )
  
  # Predictions using BEST model per curve (skip NULL preds safely)
  preds_tbl <- models_tbl %>%
    mutate(
      pred_grid = map(data, make_pred_grid, n = 120),
      preds = map2(best_model, pred_grid, augment_safe)
    ) %>%
    select(Id, Rep, SalinityLevel, Stress, preds) %>%
    filter(!map_lgl(preds, is.null)) %>%
    unnest(preds)
  
  # Curve fit plot per ID
  p_fit <- ggplot() +
    geom_point(data = df_id, aes(t, log10_od_cor, col = as.factor(Rep)), alpha = 0.8) +
    { if (nrow(preds_tbl) > 0) geom_line(data = preds_tbl, aes(t, .fitted, col = as.factor(Rep)), linewidth = 0.8) } +
    facet_grid(SalinityLevel ~ Stress, scales = "free") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(this_id, " (best model per curve)"),
      x = time_label, y = "log10(OD)"
    )
  
  ggsave(
    filename = file.path(DIR_FIG, paste0(this_id, "_CurveFits_bestModel_Salinity.pdf")),
    plot = p_fit, width = 12, height = 10
  )
  
  # Parameter table (SAFE)
  params_tbl <- models_tbl %>%
    mutate(params = map(best_model, tidy_safe)) %>%
    select(Id, Rep, SalinityLevel, Stress, best_model_name, best_aicc, params) %>%
    filter(!map_lgl(params, is.null)) %>%
    unnest(params)
  
  write.csv(
    params_tbl,
    file = file.path(DIR_TAB, paste0(this_id, "_BestModel_Params_Salinity.csv")),
    row.names = FALSE
  )
  
  best_counts <- models_tbl %>%
    count(best_model_name) %>%
    mutate(
      perc = round(100 * n / sum(n), 1),
      Id = this_id
    )
  
  all_summaries[[this_id]] <- best_counts
}

# Winner summary across IDs
if (length(all_summaries) > 0) {
  winner_summary <- bind_rows(all_summaries) %>%
    arrange(Id, desc(n))
  
  write.csv(
    winner_summary,
    file = file.path(DIR_TAB, "BestModel_WinnerSummary_allIDs_Salinity.csv"),
    row.names = FALSE
  )
}

message("Main loop done. Outputs in:\n  ", P_OUT(""))

# =========================================================
# SUPPLEMENTARY: Fig_S2.pdf (ALL IDs in ONE PDF)
# - selects low / opt / high salinity per Id x Stress (opt = max median mumax)
# - legend shows numeric salinity values
# =========================================================
message("Building Fig_S2.pdf: low/opt/high salinity per ID × Stress (numeric legend) ...")

rds_files <- file.path(DIR_RDS, paste0(ids_to_run, "_LogisticMod_Salinity.rds"))
rds_files <- rds_files[file.exists(rds_files)]

if (length(rds_files) == 0) {
  warning("No RDS model files found in:\n  ", DIR_RDS,
          "\nRun the main loop first so *_LogisticMod_Salinity.rds files exist.")
} else {
  
  all_models <- purrr::map_dfr(rds_files, readRDS)
  
  # Predictions for best models
  all_preds <- all_models %>%
    mutate(
      pred_grid = purrr::map(data, make_pred_grid, n = 120),
      preds = purrr::map2(best_model, pred_grid, augment_safe)
    ) %>%
    select(Id, Rep, SalinityLevel, Stress, preds) %>%
    filter(!map_lgl(preds, is.null)) %>%
    tidyr::unnest(preds) %>%
    mutate(
      Rep = as.factor(Rep),
      Stress = as.factor(Stress),
      SalinityLevel = as.numeric(SalinityLevel)
    )
  
  # mumax per curve -> pick opt salinity per Id x Stress
  mumax_tbl <- all_models %>%
    mutate(params = purrr::map(best_model, tidy_safe)) %>%
    select(Id, Rep, SalinityLevel, Stress, params) %>%
    filter(!map_lgl(params, is.null)) %>%
    tidyr::unnest(params) %>%
    filter(term == "mumax") %>%
    transmute(
      Id = as.character(Id),
      Rep = as.character(Rep),
      SalinityLevel = as.numeric(SalinityLevel),
      Stress = as.character(Stress),
      mumax = as.numeric(estimate)
    ) %>%
    filter(is.finite(mumax))
  
  if (nrow(mumax_tbl) == 0) {
    warning("No mumax estimates found (term == 'mumax'). Cannot pick opt salinity.")
  } else {
    
    sal_pick <- mumax_tbl %>%
      group_by(Id, Stress, SalinityLevel) %>%
      summarise(med_mumax = median(mumax, na.rm = TRUE), .groups = "drop") %>%
      group_by(Id, Stress) %>%
      summarise(
        low_sal  = min(SalinityLevel, na.rm = TRUE),
        high_sal = max(SalinityLevel, na.rm = TRUE),
        opt_sal  = SalinityLevel[which.max(med_mumax)],
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = c(low_sal, opt_sal, high_sal),
        names_to = "sal_cat",
        values_to = "SalinityLevel"
      ) %>%
      distinct(Id, Stress, SalinityLevel) %>%
      mutate(SalLabel = as.character(SalinityLevel))   # numeric legend
    
    d_plot <- d_run %>%
      mutate(
        Id = as.character(Id),
        Stress = as.character(Stress),
        Rep = as.factor(Rep),
        SalinityLevel = as.numeric(SalinityLevel)
      ) %>%
      inner_join(sal_pick, by = c("Id", "Stress", "SalinityLevel")) %>%
      mutate(SalLabel = factor(SalLabel, levels = sort(unique(SalLabel))))
    
    preds_plot <- all_preds %>%
      mutate(
        Id = as.character(Id),
        Stress = as.character(Stress),
        SalinityLevel = as.numeric(SalinityLevel)
      ) %>%
      inner_join(sal_pick, by = c("Id", "Stress", "SalinityLevel")) %>%
      mutate(SalLabel = factor(SalLabel, levels = levels(d_plot$SalLabel)))
    
    if (nrow(d_plot) == 0 || nrow(preds_plot) == 0) {
      warning("After filtering to low/opt/high salinity, nothing left to plot.")
    } else {
      
      # stable ordering
      id_levels <- ids_to_run
      stress_levels <- sort(unique(d_plot$Stress))
      
      d_plot$Id <- factor(d_plot$Id, levels = id_levels)
      preds_plot$Id <- factor(preds_plot$Id, levels = id_levels)
      
      d_plot$Stress <- factor(d_plot$Stress, levels = stress_levels)
      preds_plot$Stress <- factor(preds_plot$Stress, levels = stress_levels)
      
      p_s2 <- ggplot() +
        geom_point(
          data = d_plot,
          aes(x = t, y = log10_od_cor, colour = Rep),
          alpha = 0.45, size = 0.7
        ) +
        geom_line(
          data = preds_plot,
          aes(x = t, y = .fitted, colour = Rep, linetype = SalLabel),
          linewidth = 0.6
        ) +
        facet_grid(
          rows = vars(Id),
          cols = vars(Stress),
          scales = "free"
        ) +
        theme_bw(base_size = 10) +
        labs(
          title = "Monoculture growth curves: best-model fits (all IDs) - selected salinity levels per pH",
          x = time_label,
          y = "log10(OD)",
          linetype = "Salinity",
          colour = "Rep"
        ) +
        theme(
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
        )
      
      fig_s2_file <- file.path(DIR_FIG, "Fig_S2.pdf")
      ggsave(
        filename = fig_s2_file,
        plot = p_s2,
        width = 14,
        height = max(6, 0.9 * length(ids_to_run) + 3)
      )
      ggsave(
        filename = file.path(DIR_FIG, "Fig_S2.png"),
        plot = p_s2,
        width = 14,
        height = max(6, 0.9 * length(ids_to_run) + 3),
        dpi = 300
      )

      message("Saved: ", fig_s2_file)
    }
  }
}


### Logistic/Growth models for pH experiment (project-root clean pipeline)
### - SAME curve-fit logic as your TEMP script:
###   baranyi, baranyi_without_lag, gompertz, buchanan, buchanan_without_lag
###   best per curve via AICc
### - Uses scripts/growth_curve_models.R (your existing file)
### - Project-root paths:
###     data/...      for inputs
###     results/...   for outputs
### - Reproducible per-curve starts (seeded by Id×Rep×TempLevel×pH_Level)
### - SUPPLEMENTARY: Fig_S3.pdf (ALL IDs; low/opt/high pH per Id×TempLevel)
###   legend shows NUMERIC pH values (not low/opt/high)
### - Multi-page supplementary PDF is OFF

rm(list = ls())

suppressPackageStartupMessages({
  library(nls.multstart)
  library(ggplot2)
  library(broom)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(MuMIn)     # AICc
  library(readr)
  library(tibble)
})

# -------------------------
# 0) Project-root paths
# -------------------------
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

P_IN  <- function(...) file.path(ROOT, "data", ...)
P_SRC <- function(...) file.path(ROOT, "scripts", ...)
P_OUT <- function(...) file.path(ROOT, "results", ...)

DIR_RDS <- P_OUT("rds")
DIR_FIG <- P_OUT("figures")
DIR_TAB <- P_OUT("tables")

dir.create(DIR_RDS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 0b) Input + model source
# -------------------------
data_path <- P_IN("Cut_OD_data", "cutPoint_D17_edit2_pH.csv")
model_src <- P_SRC("06_growth_curve_models.R")

if (!file.exists(data_path)) {
  stop("data_path not found:\n  ", data_path,
       "\n\nExpected:\n  data/Cut_OD_data/cutPoint_D17_edit2_pH.csv",
       "\nCurrent getwd():\n  ", getwd())
}
if (!file.exists(model_src)) {
  stop("model_src not found:\n  ", model_src,
       "\n\nExpected:\n  scripts/growth_curve_models.R",
       "\nCurrent getwd():\n  ", getwd())
}

source(model_src)

# ---- TIME IN MINUTES ----
# If your raw column t is already minutes, keep this as 1.
# If your raw column t is hours, set time_multiplier <- 60.
time_multiplier <- 1
time_label <- "Time (min)"

# -------------------------
# 1) Load + clean data
#    Expect columns: t, OD, pH_Level, TempLevel, Rep OR Replicate, Id
# -------------------------
d_raw <- readr::read_csv(data_path, show_col_types = FALSE) %>%
  select(-matches("^X\\.[0-9]+$"), -any_of("X"))

req_base <- c("OD", "t", "pH_Level", "TempLevel", "Id")
missing_base <- setdiff(req_base, names(d_raw))
if (length(missing_base) > 0) {
  stop("Missing required columns in CSV: ", paste(missing_base, collapse = ", "))
}

rep_col <- if ("Rep" %in% names(d_raw)) "Rep" else if ("Replicate" %in% names(d_raw)) "Replicate" else NA_character_
if (is.na(rep_col)) stop("CSV must contain either 'Rep' or 'Replicate' column.")

d <- d_raw %>%
  transmute(
    t  = as.numeric(t) * time_multiplier,
    OD = as.numeric(OD),
    pH_Level  = as.numeric(pH_Level),
    TempLevel = as.numeric(TempLevel),
    Rep = as.factor(.data[[rep_col]]),
    Id  = as.character(Id),
    
    od_cor2 = pmax(OD, 0.040),
    log10_od_cor = log10(od_cor2)
  ) %>%
  filter(
    !is.na(log10_od_cor), is.finite(log10_od_cor),
    !is.na(t),            is.finite(t),
    !is.na(pH_Level),     is.finite(pH_Level),
    !is.na(TempLevel),    is.finite(TempLevel),
    !is.na(Rep),
    !is.na(Id),
    Id != "BLANK"
  ) %>%
  arrange(Id, TempLevel, pH_Level, Rep, t)

# -------------------------
# 2) Deterministic seeding per curve (Id×Rep×TempLevel×pH_Level)
# -------------------------
seed_from_group <- function(Id, Rep, TempLevel, pH_Level) {
  x <- paste(Id, Rep, TempLevel, pH_Level, sep = "|")
  as.integer(sum(utf8ToInt(x)) %% .Machine$integer.max)
}

# -------------------------
# 3) Model fitting helpers (SAME bounds / iter as your reference script)
# -------------------------
safe_nls <- function(formula, data, iter, start_lower, start_upper, lower, upper) {
  tryCatch(
    nls_multstart(
      formula = formula,
      data = data,
      iter = iter,
      start_lower = start_lower,
      start_upper = start_upper,
      supp_errors = "Y",
      na.action = na.omit,
      lower = lower,
      upper = upper
    ),
    error = function(e) NULL
  )
}

fit_all_models_one_curve <- function(df_curve, Id, Rep, TempLevel, pH_Level) {
  
  if (nrow(df_curve) < 8) {
    return(list(
      fits_baranyi = NULL,
      fits_baranyi_without_lag = NULL,
      fits_gompertz = NULL,
      fits_buchanan = NULL,
      fits_buchanan_without_lag = NULL
    ))
  }
  
  set.seed(seed_from_group(Id, Rep, TempLevel, pH_Level))
  
  list(
    fits_baranyi = safe_nls(
      log10_od_cor ~ baranyi(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_baranyi_without_lag = safe_nls(
      log10_od_cor ~ baranyi_without_lag(log10_nmax, log10_n0, mumax, t = t),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10)
    ),
    
    fits_gompertz = safe_nls(
      log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_buchanan = safe_nls(
      log10_od_cor ~ buchanan(log10_nmax, log10_n0, mumax, t = t, lag),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0, lag = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5, lag = 1500),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0, lag = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10, lag = 3000)
    ),
    
    fits_buchanan_without_lag = safe_nls(
      log10_od_cor ~ buchanan_without_lag(log10_nmax, log10_n0, mumax, t = t),
      data = df_curve, iter = 500,
      start_lower = c(log10_nmax = -3,  log10_n0 = -3, mumax = 0),
      start_upper = c(log10_nmax = -0.2,log10_n0 = -1, mumax = 5),
      lower =       c(log10_nmax = -5,  log10_n0 = -5, mumax = 0),
      upper =       c(log10_nmax =  0,  log10_n0 =  0, mumax = 10)
    )
  )
}

aicc_or_na <- function(model_obj) {
  if (is.null(model_obj)) return(NA_real_)
  tryCatch(MuMIn::AICc(model_obj), error = function(e) NA_real_)
}

pick_best_model <- function(model_list) {
  aiccs <- purrr::map_dbl(model_list, aicc_or_na)
  if (all(is.na(aiccs))) {
    return(list(best_name = NA_character_, best_model = NULL, best_aicc = NA_real_))
  }
  best_name <- names(which.min(aiccs))
  list(
    best_name = best_name,
    best_model = model_list[[best_name]],
    best_aicc = min(aiccs, na.rm = TRUE)
  )
}

make_pred_grid <- function(df, n = 120) {
  data.frame(t = seq(min(df$t, na.rm = TRUE), max(df$t, na.rm = TRUE), length.out = n))
}

augment_safe <- function(model_obj, newdata) {
  if (is.null(model_obj)) return(NULL)
  tryCatch(broom::augment(model_obj, newdata = newdata), error = function(e) NULL)
}

# SAFE tidy (avoids chol2inv crash)
tidy_safe <- function(model_obj) {
  if (is.null(model_obj)) return(NULL)
  
  td <- tryCatch(broom::tidy(model_obj), error = function(e) NULL)
  if (!is.null(td)) return(td)
  
  cf <- tryCatch(coef(model_obj), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  
  tibble(
    term = names(cf),
    estimate = as.numeric(cf),
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_
  )
}

# -------------------------
# 4) Run per-ID (multi-page supplementary PDF OFF)
# -------------------------
ids_to_run <- unique(d$Id)

exclude_rules <- tibble::tribble(
  ~Id, ~TempLevel, ~pH_Level, ~Rep
  # leave empty unless needed
)

apply_exclusions <- function(df, rules) {
  if (nrow(rules) == 0) return(df)
  df %>%
    mutate(Rep = as.character(Rep)) %>%
    anti_join(
      rules %>% mutate(Rep = as.character(Rep)),
      by = c("Id", "TempLevel", "pH_Level", "Rep")
    )
}

d_run <- d %>% apply_exclusions(exclude_rules)

all_summaries <- list()

for (this_id in ids_to_run) {
  
  message("----- Running ID: ", this_id, " -----")
  
  df_id <- d_run %>% filter(Id == this_id)
  
  if (nrow(df_id) == 0) {
    message("No data for ", this_id, " (skipping).")
    next
  }
  
  # Fit per curve: Id × Rep × TempLevel × pH_Level
  models_tbl <- df_id %>%
    group_by(Id, Rep, TempLevel, pH_Level) %>%
    nest() %>%
    mutate(
      fits = pmap(list(data, Id, Rep, TempLevel, pH_Level),
                  ~ fit_all_models_one_curve(..1, ..2, ..3, ..4, ..5)),
      best = map(fits, pick_best_model),
      best_model_name = map_chr(best, "best_name"),
      best_aicc = map_dbl(best, "best_aicc"),
      best_model = map(best, "best_model")
    ) %>%
    ungroup()
  
  rds_file <- file.path(DIR_RDS, paste0(this_id, "_LogisticMod_pH.rds"))
  saveRDS(models_tbl, rds_file)
  message("Saved: ", rds_file)
  
  # AICc histogram
  model_stack <- models_tbl %>%
    select(Id, Rep, TempLevel, pH_Level, fits) %>%
    mutate(
      fits_baranyi = map(fits, "fits_baranyi"),
      fits_baranyi_without_lag = map(fits, "fits_baranyi_without_lag"),
      fits_buchanan = map(fits, "fits_buchanan"),
      fits_buchanan_without_lag = map(fits, "fits_buchanan_without_lag"),
      fits_gompertz = map(fits, "fits_gompertz")
    ) %>%
    select(-fits) %>%
    pivot_longer(
      cols = starts_with("fits_"),
      names_to = "model",
      values_to = "model_obj"
    ) %>%
    mutate(
      aic = map_dbl(model_obj, aicc_or_na),
      model_name = case_when(
        model == "fits_baranyi" ~ "(a) baranyi",
        model == "fits_baranyi_without_lag" ~ "(b) baranyi without lag",
        model == "fits_buchanan" ~ "(c) buchanan",
        model == "fits_buchanan_without_lag" ~ "(d) buchanan without lag",
        model == "fits_gompertz" ~ "(e) gompertz",
        TRUE ~ model
      )
    )
  
  aic_means <- model_stack %>%
    filter(!is.na(aic)) %>%
    group_by(model, model_name) %>%
    summarise(
      mean_aic = mean(aic),
      median_aic = median(aic),
      n = n(),
      .groups = "drop"
    )
  
  p_aic <- ggplot(model_stack, aes(aic)) +
    geom_histogram(fill = "lightgrey", col = "black", binwidth = 10) +
    geom_vline(data = aic_means, aes(xintercept = mean_aic), col = "red") +
    geom_vline(data = aic_means, aes(xintercept = median_aic), col = "blue") +
    facet_wrap(~ model_name, ncol = 2) +
    theme_bw(base_size = 14) +
    ylab("Count") + xlab("AICc score") +
    theme(strip.background = element_blank(),
          strip.text = element_text(hjust = 0))
  
  ggsave(
    filename = file.path(DIR_FIG, paste0(this_id, "_AICc_hist_pH.pdf")),
    plot = p_aic, width = 10, height = 8
  )
  
  # Predictions using BEST model per curve
  preds_tbl <- models_tbl %>%
    mutate(
      pred_grid = map(data, make_pred_grid, n = 120),
      preds = map2(best_model, pred_grid, augment_safe)
    ) %>%
    select(Id, Rep, TempLevel, pH_Level, preds) %>%
    filter(!map_lgl(preds, is.null)) %>%
    unnest(preds)
  
  # Curve plot per ID
  p_fit <- ggplot() +
    geom_point(data = df_id, aes(t, log10_od_cor, col = as.factor(Rep)), alpha = 0.8) +
    { if (nrow(preds_tbl) > 0) geom_line(data = preds_tbl, aes(t, .fitted, col = as.factor(Rep)), linewidth = 0.8) } +
    facet_grid(TempLevel ~ pH_Level, scales = "free") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(this_id, " (best model per curve)"),
      x = time_label, y = "log10(OD)"
    )
  
  ggsave(
    filename = file.path(DIR_FIG, paste0(this_id, "_CurveFits_bestModel_pH.pdf")),
    plot = p_fit, width = 12, height = 10
  )
  
  # Parameter table (SAFE)
  params_tbl <- models_tbl %>%
    mutate(params = map(best_model, tidy_safe)) %>%
    select(Id, Rep, TempLevel, pH_Level, best_model_name, best_aicc, params) %>%
    filter(!map_lgl(params, is.null)) %>%
    unnest(params)
  
  write.csv(
    params_tbl,
    file = file.path(DIR_TAB, paste0(this_id, "_BestModel_Params_pH.csv")),
    row.names = FALSE
  )
  
  best_counts <- models_tbl %>%
    count(best_model_name) %>%
    mutate(
      perc = round(100 * n / sum(n), 1),
      Id = this_id
    )
  
  all_summaries[[this_id]] <- best_counts
}

# Winner summary across IDs
if (length(all_summaries) > 0) {
  winner_summary <- bind_rows(all_summaries) %>%
    arrange(Id, desc(n))
  
  write.csv(
    winner_summary,
    file = file.path(DIR_TAB, "BestModel_WinnerSummary_allIDs_pH.csv"),
    row.names = FALSE
  )
}

message("Main loop done. Outputs in:\n  ", P_OUT(""))

# =========================================================
# SUPPLEMENTARY: Fig_S3.pdf (ALL IDs in ONE PDF)
# - low / opt / high pH per Id x TempLevel (opt = max median mumax)
# - legend shows numeric pH values
# =========================================================
message("Building Fig_S3.pdf: low/opt/high pH per ID × TempLevel (numeric legend) ...")

rds_files <- file.path(DIR_RDS, paste0(ids_to_run, "_LogisticMod_pH.rds"))
rds_files <- rds_files[file.exists(rds_files)]

if (length(rds_files) == 0) {
  warning("No RDS model files found in:\n  ", DIR_RDS,
          "\nRun the main loop first so *_LogisticMod_pH.rds files exist.")
} else {
  
  all_models <- purrr::map_dfr(rds_files, readRDS)
  
  # predictions for best models
  all_preds <- all_models %>%
    mutate(
      pred_grid = purrr::map(data, make_pred_grid, n = 120),
      preds = purrr::map2(best_model, pred_grid, augment_safe)
    ) %>%
    select(Id, Rep, TempLevel, pH_Level, preds) %>%
    filter(!map_lgl(preds, is.null)) %>%
    tidyr::unnest(preds) %>%
    mutate(
      Rep = as.factor(Rep),
      TempLevel = as.numeric(TempLevel),
      pH_Level = as.numeric(pH_Level)
    )
  
  # mumax per curve -> pick opt pH per Id x TempLevel
  mumax_tbl <- all_models %>%
    mutate(params = purrr::map(best_model, tidy_safe)) %>%
    select(Id, Rep, TempLevel, pH_Level, params) %>%
    filter(!map_lgl(params, is.null)) %>%
    tidyr::unnest(params) %>%
    filter(term == "mumax") %>%
    transmute(
      Id = as.character(Id),
      Rep = as.character(Rep),
      TempLevel = as.numeric(TempLevel),
      pH_Level = as.numeric(pH_Level),
      mumax = as.numeric(estimate)
    ) %>%
    filter(is.finite(mumax))
  
  if (nrow(mumax_tbl) == 0) {
    warning("No mumax estimates found (term == 'mumax'). Cannot pick opt pH.")
  } else {
    
    ph_pick <- mumax_tbl %>%
      group_by(Id, TempLevel, pH_Level) %>%
      summarise(med_mumax = median(mumax, na.rm = TRUE), .groups = "drop") %>%
      group_by(Id, TempLevel) %>%
      summarise(
        low_pH  = min(pH_Level, na.rm = TRUE),
        high_pH = max(pH_Level, na.rm = TRUE),
        opt_pH  = pH_Level[which.max(med_mumax)],
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = c(low_pH, opt_pH, high_pH),
        names_to = "ph_cat",
        values_to = "pH_Level"
      ) %>%
      distinct(Id, TempLevel, pH_Level) %>%
      mutate(pH_label = sprintf("%.1f", pH_Level))  # numeric legend
    
    d_plot <- d_run %>%
      mutate(
        Id = as.character(Id),
        TempLevel = as.numeric(TempLevel),
        Rep = as.factor(Rep),
        pH_Level = as.numeric(pH_Level)
      ) %>%
      inner_join(ph_pick, by = c("Id", "TempLevel", "pH_Level")) %>%
      mutate(pH_label = factor(pH_label, levels = sort(unique(pH_label))))
    
    preds_plot <- all_preds %>%
      mutate(
        Id = as.character(Id),
        TempLevel = as.numeric(TempLevel),
        pH_Level = as.numeric(pH_Level)
      ) %>%
      inner_join(ph_pick, by = c("Id", "TempLevel", "pH_Level")) %>%
      mutate(pH_label = factor(pH_label, levels = levels(d_plot$pH_label)))
    
    if (nrow(d_plot) == 0 || nrow(preds_plot) == 0) {
      warning("After filtering to low/opt/high pH, nothing left to plot.")
    } else {
      
      # stable facet ordering
      id_levels <- ids_to_run
      temp_levels <- sort(unique(d_plot$TempLevel))
      
      d_plot$Id <- factor(d_plot$Id, levels = id_levels)
      preds_plot$Id <- factor(preds_plot$Id, levels = id_levels)
      
      d_plot$TempLevel <- factor(d_plot$TempLevel, levels = temp_levels)
      preds_plot$TempLevel <- factor(preds_plot$TempLevel, levels = temp_levels)
      
      p_s3 <- ggplot() +
        geom_point(
          data = d_plot,
          aes(x = t, y = log10_od_cor, colour = Rep),
          alpha = 0.45, size = 0.7
        ) +
        geom_line(
          data = preds_plot,
          aes(x = t, y = .fitted, colour = Rep, linetype = pH_label),
          linewidth = 0.6
        ) +
        facet_grid(
          rows = vars(Id),
          cols = vars(TempLevel),
          scales = "free"
        ) +
        theme_bw(base_size = 10) +
        labs(
          title = "Growth curves: best-model fits (all IDs) - selected pH levels per temperature",
          x = time_label,
          y = "log10(OD)",
          linetype = "pH",
          colour = "Rep"
        ) +
        theme(
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
        )
      
      fig_s3_file <- file.path(DIR_FIG, "Fig_S3.pdf")
      ggsave(
        filename = fig_s3_file,
        plot = p_s3,
        width = 14,
        height = max(6, 0.9 * length(ids_to_run) + 3)
      )
      ggsave(
        filename = file.path(DIR_FIG, "Fig_S3.png"),
        plot = p_s3,
        width = 14,
        height = max(6, 0.9 * length(ids_to_run) + 3),
        dpi = 300
      )

      message("Saved: ", fig_s3_file)
    }
  }
}
