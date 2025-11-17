# ======================================================================
# Fig 1: Diversity, Community Composition, and Environmental Gradients
#
# Products:
# - Fig 1a: Shannon diversity (boxplots + significance)
# - Fig 1b: db-RDA of community composition vs pH, salinity, temperature
# - Fig 1c: PERMANOVA pairwise R² heatmap
# - Table S4: Tukey post-hoc contrasts + effect sizes for Shannon diversity
# - Table S5: Pairwise PERMANOVA on Bray–Curtis distances
#
# Overview:
# - Loads preprocessed ASV abundance matrix and aligned metadata
# - Recodes stress treatments and builds Bray–Curtis distance matrix
# - Performs db-RDA of community composition along pH, salinity, and temperature
# - Runs global PERMANOVA and pairwise PERMANOVA by Stress
# - Computes Shannon diversity, runs ANOVA + post-hoc tests, and effect sizes
# - Generates and saves all corresponding plots and tables for Fig 1 & supplements
#
# Folder conventions:
# - Inputs (cleaned):  rds/       (via P_RDS())
# - Plots:             plots/     (via P_FIG())
# - Tables (CSV):      tables/    (via P_TAB())
#
# Outline:
# 1) Load cleaned ASV matrix and metadata
# 2) Recode Stress factor and compute Bray–Curtis distances
# 3) db-RDA on pH × Salinity × Temperature and save Fig 1b
# 4) Global + pairwise PERMANOVA; build R² heatmap and save Fig 1c & Table S5
# 5) Compute Shannon diversity, run ANOVA + post-hoc & effect sizes
#    and save Fig 1a, Table S4, and Shannon_ANOVA_effectsizes.csv

# ======================================================================
# Fig 1 (diversity, PERMANOVA, dbRDA) -> Tables S4–S5, Fig 1a/b/c
source("scripts/utils_functions.R"); ensure_packages()

ASV_t_clean        <- readRDS(P_RDS("ASV_matrix_clean.rds"))
env_option_A_clean <- readRDS(P_RDS("metadata_clean.rds"))

# Stress factor & palette
env_option_A_clean$Stress <- env_option_A_clean$Stress %>% as.character() %>%
  {replace(., . == "Salinity", "Sal")} %>% {replace(., . == "Temperature", "Temp")} %>%
  factor(levels = levels_short)

# Distance matrix
Distances <- vegan::vegdist(ASV_t_clean, method = "bray")

# --- db-RDA (Fig 1b) ---
env_dbRDA <- env_option_A_clean %>%
  mutate(pH_num = as.numeric(as.character(pH)),
         Salinity_num = as.numeric(as.character(Salinity)),
         Temp_num = as.numeric(as.character(Temp)))

dbRDA_grad <- vegan::capscale(Distances ~ pH_num * Salinity_num * Temp_num, data = env_dbRDA)
site_scores_grad <- vegan::scores(dbRDA_grad, display = "sites") %>% as.data.frame()
site_scores_grad$Stress <- env_option_A_clean$Stress
var_expl_grad <- { eig <- vegan::eigenvals(dbRDA_grad, constrained = TRUE); round(100 * eig / sum(eig), 1) }

dbRDA_plot_grad <- ggplot(site_scores_grad, aes(CAP1, CAP2, colour = Stress)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Stress, colour = Stress), type = "t", level = 0.95, linewidth = 0.8, alpha = 0.25) +
  labs(x = paste0("RDA1 (", var_expl_grad[1], "%)"), y = paste0("RDA2 (", var_expl_grad[2], "%)")) +
  scale_colour_manual(values = pal_short, limits = levels_short, breaks = levels_short) +
  theme_classic(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14),
        legend.position = "right")

ggsave(P_FIG("Fig_1b_dbRDA_gradients.tiff"), label_plot(dbRDA_plot_grad, "1b"),
       width = 18, height = 15, units = "cm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")

# --- PERMANOVA + pairwise heatmap (Fig 1c & Table S5) ---
global_perm <- vegan::adonis2(Distances ~ Stress, data = env_option_A_clean, permutations = 999)
cat("Global PERMANOVA: R² =", round(global_perm$R2[1], 2), ", p =", signif(global_perm$`Pr(>F)`[1], 3), "\n")

pw <- pairwiseAdonis::pairwise.adonis(x = Distances, factors = env_option_A_clean$Stress, perm = 999) %>% as.data.frame()
pw$Group1 <- sub(" vs .*", "", pw$pairs); pw$Group2 <- sub(".* vs ", "", pw$pairs)
present_short <- levels_short[levels_short %in% unique(as.character(env_option_A_clean$Stress))]
mat <- matrix(NA_real_, nrow = length(present_short), ncol = length(present_short),
              dimnames = list(present_short, present_short))
for (i in seq_len(nrow(pw))) { g1 <- pw$Group1[i]; g2 <- pw$Group2[i]; if (g1 %in% present_short && g2 %in% present_short) mat[g1, g2] <- pw$R2[i] }
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
mat_long <- reshape2::melt(mat, varnames = c("Group1","Group2"), value.name = "R2", na.rm = TRUE)

pcol <- intersect(c("p.adjusted","p.adjustment","p.adj"), names(pw))
if (length(pcol) == 0) stop("Could not find adjusted p-value column in pairwise.adonis output.")
sig_map <- pw[, c("Group1","Group2", pcol[1])]; names(sig_map)[3] <- "p.adj"
sig_map$stars <- cut(sig_map$p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, 1), labels = c("***","**","*",""))
mat_long <- mat_long %>% left_join(sig_map, by = c("Group1","Group2"))
pretty_present <- unname(pretty_map[present_short])
mat_long <- mat_long %>% mutate(
  Group1 = dplyr::recode(Group1, !!!as.list(pretty_map)),
  Group2 = dplyr::recode(Group2, !!!as.list(pretty_map)),
  Group1 = factor(Group1, levels = pretty_present),
  Group2 = factor(Group2, levels = pretty_present)
)

heatmap_plot <- ggplot(mat_long, aes(x = Group1, y = Group2, fill = R2)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkred", na.value = "grey90") +
  geom_text(aes(label = stars), vjust = 0.5, size = 5, color = "black", na.rm = TRUE) +
  labs(x = NULL, y = NULL, fill = expression(R^2)) +
  theme_classic(base_family = "Helvetica") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14))

ggsave(P_FIG("Fig_1c_PERMANOVA_pairwise_heatmap.tiff"), label_plot(heatmap_plot, "1c"),
       width = 18, height = 15, units = "cm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")

readr::write_csv(pw %>% mutate(across(where(is.numeric), ~round(., 3))), P_TAB("Table_S5_PERMANOVA_pairwise.csv"))

# --- Shannon diversity (Fig 1a & Table S4 + Shannon_ANOVA) ---
if (!"Id" %in% names(env_option_A_clean)) { env_option_A_clean <- env_option_A_clean %>% tibble::rownames_to_column("Id") }

grab_col <- function(df, candidates, default = NA_real_) { for (nm in candidates) if (nm %in% names(df)) return(df[[nm]]); rep(default, nrow(df)) }

Data_richness <- vegan::diversity(ASV_t_clean, index = "shannon")
Data_richness_df <- tibble::tibble(Id = names(Data_richness), shan_div = as.numeric(Data_richness))

div_df <- env_option_A_clean %>%
  dplyr::select(Id, Stress) %>%
  dplyr::left_join(Data_richness_df, by = "Id") %>%
  dplyr::mutate(
    Stress        = factor(as.character(Stress), levels = levels_short),
    Stress_pretty = dplyr::recode(as.character(Stress), !!!pretty_map, .default = as.character(Stress)),
    Stress_pretty = factor(Stress_pretty, levels = pretty_levels)
  )

fit_shan      <- stats::lm(shan_div ~ Stress, data = div_df)
Shannon_anova <- stats::aov(shan_div ~ Stress, data = div_df)

tukey_df <- TukeyHSD(Shannon_anova, conf.level = 0.95)[[1]] %>% as.data.frame() %>% tibble::rownames_to_column("pairwise")
names(tukey_df) <- make.names(names(tukey_df))

plot_shandiv_box <- ggplot(div_df, aes(x = Stress_pretty, y = shan_div, fill = Stress_pretty)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.6) +
  scale_fill_manual(values = pal_pretty, limits = pretty_levels) +
  scale_x_discrete(labels = scales::label_wrap(12)) +
  labs(x = NULL, y = "Shannon diversity (H′)") +
  theme_classic(base_family = "Helvetica") +
  theme(axis.text.x  = element_text(angle = 50, hjust = 1, vjust = 1, size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        plot.margin = margin(8, 12, 8, 8))

assump_shapiro <- shapiro.test(residuals(fit_shan))
assump_levene  <- car::leveneTest(shan_div ~ Stress, div_df)
welch_alt      <- oneway.test(shan_div ~ Stress, div_df, var.equal = FALSE)
anova_type2 <- car::Anova(fit_shan, type = 2) %>% broom::tidy()

eta2_raw  <- effectsize::eta_squared(fit_shan, partial = FALSE, ci = 0.95) %>% as.data.frame()
omega_raw <- effectsize::omega_squared(fit_shan, partial = FALSE) %>% as.data.frame()

eta2_tbl <- tibble::tibble(
  term   = grab_col(eta2_raw, c("Effect","Term","Parameter"), default = NA_character_),
  eta2   = grab_col(eta2_raw, c("Eta2_partial","Eta2")),
  ci_low = grab_col(eta2_raw, c("CI_low","CI_low_partial","CI_low_Eta2")),
  ci_high= grab_col(eta2_raw, c("CI_high","CI_high_partial","CI_high_Eta2"))
) %>% dplyr::distinct()

omega2_tbl <- tibble::tibble(
  term   = grab_col(omega_raw, c("Effect","Term","Parameter"), default = NA_character_),
  omega2 = grab_col(omega_raw, c("Omega2_partial","Omega2"))
) %>% dplyr::distinct()

S1a <- anova_type2 %>% dplyr::filter(term == "Stress") %>%
  dplyr::transmute(term, df = as.integer(df), ss = sumsq, ms = sumsq/df,
                   F = statistic, p = p.value, N = nrow(div_df), k = dplyr::n_distinct(div_df$Stress)) %>%
  dplyr::left_join(eta2_tbl,  by = "term") %>%
  dplyr::left_join(omega2_tbl, by = "term") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))
readr::write_csv(S1a, P_TAB("Shannon_ANOVA_effectsizes.csv"))

emm <- emmeans::emmeans(fit_shan, ~ Stress)
pw_contr <- emmeans::contrast(emm, method = "pairwise")
tuk_raw  <- summary(pw_contr, infer = c(FALSE, TRUE), adjust = "tukey") %>% as.data.frame()
stat_col <- if ("t.ratio" %in% names(tuk_raw)) "t.ratio" else if ("z.ratio" %in% names(tuk_raw)) "z.ratio" else NA_character_
tuk_tbl <- tuk_raw %>% dplyr::transmute(
  contrast = .data$contrast, estimate = .data$estimate, SE = .data$SE,
  df = if ("df" %in% names(tuk_raw)) .data$df else NA_real_,
  stat = if (!is.na(stat_col)) .data[[stat_col]] else NA_real_,
  p.adj = if ("p.value" %in% names(tuk_raw)) .data$p.value else NA_real_
)

sigma_hat <- sigma(fit_shan); df_resid  <- df.residual(fit_shan); J_corr <- 1 - 3/(4*df_resid - 1)
S1b <- tuk_tbl %>% dplyr::mutate(
  hedges_g = J_corr * (estimate / sigma_hat),
  SE_g     = J_corr * (SE / sigma_hat),
  t_crit   = ifelse(is.finite(df), qt(0.975, df), NA_real_),
  g_low    = hedges_g - t_crit * SE_g,
  g_high   = hedges_g + t_crit * SE_g
) %>% dplyr::select(contrast, estimate, SE, df, stat, p.adj, hedges_g, g_low, g_high) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3))) %>% dplyr::arrange(p.adj)
readr::write_csv(S1b, P_TAB("Table_S4_Tukey_posthoc_effectsizes.csv"))

tukey_pairs <- TukeyHSD(Shannon_anova, conf.level = 0.95)[[1]] %>% as.data.frame() %>% tibble::rownames_to_column("pairwise")
names(tukey_pairs) <- make.names(names(tukey_pairs))
sig_df <- tukey_pairs %>% tidyr::separate(pairwise, into = c("g1","g2"), sep = "-") %>%
  dplyr::mutate(g1 = dplyr::recode(g1, !!!pretty_map, .default = g1),
                g2 = dplyr::recode(g2, !!!pretty_map, .default = g2)) %>%
  dplyr::filter(g1 == "Control" | g2 == "Control") %>% dplyr::arrange(p.adj) %>% dplyr::filter(p.adj < 0.05)

plot_shandiv_box_sig <- plot_shandiv_box
if (nrow(sig_df) > 0) {
  y_range <- range(div_df$shan_div, na.rm = TRUE); step <- 0.06 * diff(y_range); base <- max(div_df$shan_div, na.rm = TRUE) + step
  sig_df <- sig_df %>% dplyr::transmute(
    group1 = ifelse(g1 == "Control", g1, g2), group2 = ifelse(g1 == "Control", g2, g1),
    p.adj, y.position = base + dplyr::row_number() * step,
    p.signif = dplyr::case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", TRUE ~ "*")
  )
  plot_shandiv_box_sig <- plot_shandiv_box +
    ggpubr::stat_pvalue_manual(sig_df, label = "p.signif", xmin = "group1", xmax = "group2",
                               y.position = "y.position", tip.length = 0.01, bracket.size = 0.3, size = 3) +
    expand_limits(y = max(sig_df$y.position) + step)
}
ggsave(P_FIG("Fig_1a_Shannon_diversity_boxplot_with_signif.tiff"), label_plot(plot_shandiv_box_sig, "1a"),
       width = 18, height = 15, units = "cm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")
