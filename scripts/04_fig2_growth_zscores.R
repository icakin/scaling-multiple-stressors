# ======================================================================
# Fig 2: Growth Rate z-Scores Across Stress Combinations
#
# Products:
# - Fig 2: Compact dotplot of z-scored growth rates by taxon and stress
# - Table S7: Long-form growth rate dataset (all sources × stressors)
# - Table S8: Summary of z-scored growth (mean ± SE) per taxon × stress
#
# Overview:
# - Loads Gompertz growth-rate fits from separate Temp, pH, and Salinity runs
# - Harmonises column names and core fields (Id, growth rate, replicate, conditions)
# - Normalises raw temperature, salinity, and pH levels to experimental settings
# - Derives the Stress factor (Control, Temp, pH, Sal, and all combinations)
# - Combines all sources into a single long table and filters non-negative growth rates
# - Computes within-stress z-scores for each taxon’s growth rate
# - Summarises z-scores as mean ± SE per taxon × stress
# - Generates a multi-panel dotplot of individual and mean z-scores (Fig 2)
#
# Folder conventions:
# - Inputs (growth fits): data/    (via P_IN())
# - Tables (CSV):         tables/  (via P_TAB())
# - RDS objects:          rds/     (via P_RDS())
# - Plots:                plots/   (via P_FIG())
#
# Outline:
# 1) Load and sanitise growth-fit CSVs for Temp, pH, and Salinity experiments
# 2) Standardise core columns (taxon Id, growth rate g) and derive Stress
# 3) Bind all sources into Growth_all_sources_stressors_long and save (Table S7 + RDS)
# 4) Map taxon IDs to short labels (Genus abbreviations) for plotting
# 5) Compute z-scores of growth within each Stress and summarise (Table S8 + RDS)
# 6) Build compact dotplot (per-taxon z-scores + mean ± SE) and save Fig 2

# ======================================================================
#  Fig 2 (growth z-scores) -> Tables S7–S8, Fig 2
source("scripts/utils_functions.R"); ensure_packages()

IN_TEMP <- P_IN("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv")
IN_pH   <- P_IN("LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv")  # required
IN_SAL  <- P_IN("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv")

sanitize_names <- function(df) {
  nm <- names(df); nm <- gsub("\\p{Zs}+", " ", nm, perl = TRUE); nm <- trimws(nm)
  for (i in seq_along(nm)) if (is.na(nm[i]) || nm[i] == "") nm[i] <- paste0("V", i)
  names(df) <- make.unique(nm, sep = "_"); df
}
load_csv <- function(path){ if (!file.exists(path)) stop("Missing input: ", path); read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) |> sanitize_names() }
has_col  <- function(d, nm) nm %in% names(d)

standardise_core_cols <- function(d){
  id_col <- intersect(c("Id","OTU","Taxon","Strain","taxon","id"), names(d))[1]
  g_col  <- intersect(c("estimate","g","growth","rate","r","Growth","Rate"), names(d))[1]
  if (is.na(id_col)) stop("Could not find a taxon column among: Id/OTU/Taxon/Strain.")
  if (is.na(g_col))  stop("Could not find a growth column among: estimate/g/rate/r.")
  d %>% mutate(Id = trimf(.data[[id_col]]), g  = safe_num(.data[[g_col]]))
}

norm_temp <- function(x){ x <- safe_num(x); x[!(x %in% c(20, 38))] <- NA_real_; x }
norm_sal  <- function(x){ x <- safe_num(x); x[!(x %in% c(0, 20))]  <- NA_real_; x }
norm_pH   <- function(x){ x <- safe_num(x)
  dplyr::case_when(is.na(x) ~ NA_real_, abs(x - 5.5) <= 0.2 ~ 5.5, abs(x - 7.2) <= 0.3 ~ 7.2, TRUE ~ NA_real_) }

derive_stress <- function(temp, sal, ph){
  temp <- norm_temp(temp); sal <- norm_sal(sal); ph <- norm_pH(ph)
  ifelse(is.na(temp) | is.na(sal) | is.na(ph), NA_character_,
    { t_on  <- ifelse(temp == 38, 1L, 0L); s_on  <- ifelse(sal  == 20, 1L, 0L); p_on  <- ifelse(ph   == 5.5, 1L, 0L)
      key <- paste(p_on, s_on, t_on, sep = "")
      dplyr::recode(key,
        "000" = "Control","001" = "Temp","010" = "Sal","011" = "SalTemp",
        "100" = "pH","101" = "pHTemp","110" = "pHSal","111" = "pHSalTemp",
        .default = NA_character_) })
}

clean_one <- function(path, source_tag){
  raw <- load_csv(path) |> standardise_core_cols()
  raw %>% mutate(Source = source_tag,
    TempLevel_raw  = if (has_col(., "TempLevel"))     safe_num(TempLevel)     else NA_real_,
    Salinity_raw   = if (has_col(., "SalinityLevel")) safe_num(SalinityLevel) else NA_real_,
    pH_raw         = if (has_col(., "pH_Level"))      safe_num(pH_Level)      else NA_real_,
    TempLevel      = norm_temp(TempLevel_raw), SalinityLevel  = norm_sal(Salinity_raw), pH_Level = norm_pH(pH_raw),
    Stress         = derive_stress(TempLevel, SalinityLevel, pH_Level)) %>%
    dplyr::filter(!is.na(Stress)) %>%
    dplyr::mutate(Rep = if (has_col(., "Rep")) as.integer(Rep) else as.integer(NA)) %>%
    dplyr::select(Source, Id, g, Rep, TempLevel, SalinityLevel, pH_Level, Stress)
}

if (!file.exists(IN_pH)) stop("Required input missing: ", IN_pH)
parts <- list()
parts[["pH"]] <- clean_one(IN_pH, "pH")
if (file.exists(IN_TEMP)) parts[["Temp"]] <- clean_one(IN_TEMP, "Temp")
if (file.exists(IN_SAL))  parts[["Sal"]]  <- clean_one(IN_SAL,  "Sal")
if (length(parts) == 0L) stop("No inputs available to build the long table.")

Growth_all_sources_stressors_long <- dplyr::bind_rows(parts) %>%
  dplyr::mutate(Id = trimf(Id), Stress = factor(Stress, levels = c("Control","Temp","pH","Sal","pHTemp","SalTemp","pHSal","pHSalTemp"))) %>%
  dplyr::filter(!is.na(g), g >= 0)

Growth_all_sources_stressors_long_out <- Growth_all_sources_stressors_long %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))
readr::write_csv(Growth_all_sources_stressors_long_out, P_TAB("Table_S7_Growth_all_sources_stressors_long.csv"))
saveRDS(Growth_all_sources_stressors_long, P_RDS("Growth_all_sources_stressors_long.rds"))

taxon_map <- tibble::tribble(
  ~Id,  ~id_lab,
  "D11","Ped.", "D14","Mic.", "D17","Cur.", "I2","Yer.",
  "I8","Pse.F.", "I9","Chr.", "I11","Pse.A.", "I15","Ser.",
  "I18","Aer.", "I20","Chry.", "I22","Erw.", "I23","Jan."
)

df_lab <- Growth_all_sources_stressors_long %>% dplyr::left_join(taxon_map, by = "Id") %>% dplyr::mutate(id_lab = if_else(is.na(id_lab), Id, id_lab))

z_df <- df_lab %>% dplyr::group_by(Stress) %>% dplyr::mutate(g_z = as.numeric(scale(g))) %>% dplyr::ungroup()

Zscores_summary <- z_df %>% dplyr::group_by(Stress, id_lab) %>%
  dplyr::summarise(mu = mean(g_z, na.rm = TRUE), se = sd(g_z, na.rm = TRUE) / sqrt(sum(is.finite(g_z))), .groups = "drop")

Zscores_summary_out <- Zscores_summary %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))
readr::write_csv(Zscores_summary_out, P_TAB("Table_S8_Zscores_summary.csv"))
saveRDS(list(z_df = z_df, Zscores_summary = Zscores_summary), P_RDS("Zscores_objects.rds"))

lim_global <- max(abs(Zscores_summary$mu), na.rm = TRUE)

DOTPLOT_z_by_stress_compact <- ggplot() +
  geom_point(data = z_df %>% dplyr::left_join(Zscores_summary, by = c("Stress","id_lab")),
             aes(x = reorder_within(id_lab, mu, Stress), y = g_z),
             size = 1.1, alpha = 0.35, color = "grey30",
             position = position_jitter(width = 0.08, height = 0)) +
  geom_errorbar(data = Zscores_summary, aes(x = reorder_within(id_lab, mu, Stress), ymin = mu - se, ymax = mu + se),
                width = 0.18, linewidth = 0.3, color = "black") +
  geom_point(data = Zscores_summary, aes(x = reorder_within(id_lab, mu, Stress), y = mu, fill = mu),
             shape = 21, size = 2.8, stroke = 0.3, color = "black") +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey40") +
  facet_wrap(~ Stress, scales = "free", ncol = 4) +
  coord_flip(clip = "off") + scale_x_reordered() +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-lim_global, lim_global), oob = scales::squish, name = "z-score") +
  labs(x = NULL, y = "z-scored growth rates") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88"),
        axis.text.y = element_text(face = "italic", size = 8, margin = margin(r = 2)),
        axis.text.x = element_text(size = 8), strip.text = element_text(size = 9, face = "bold"),
        legend.position = "top", legend.direction = "horizontal", panel.spacing.x = unit(4, "pt"),
        panel.spacing.y = unit(4, "pt"), plot.margin = margin(6, 8, 6, 6))

ggsave(P_FIG("Fig_2_DOTPLOT_z_by_stress_compact.tiff"), label_plot(DOTPLOT_z_by_stress_compact, "2"),
       width = 28, height = 16, units = "cm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")
