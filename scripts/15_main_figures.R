# ======================================================================
# 15_main_figures.R
#
# Builds the two main-text figures added for the Nature Communications
# version of the manuscript. Reads only existing outputs in results/
# (no refitting, no new statistics) and writes Fig_4.png and Fig_5.png to
# results/figures/.
#
#   Fig_4.png  current Fig. 3 (panels a, b) plus the abundance
#              decomposition ladder (c; Tables S11, S19, S25) and the
#              winner-takes-all contrast (d; Table S25)
#   Fig_5.png  transfer: held-out composition error (a) and held-out
#              abundance R2 (b) across the three transfer tests
#              (Tables S15, S16, S20, S28), and leave-one-taxon-out
#              rank agreement (c) and error (d) (Table S27)
#
# Run from the repository root:  Rscript scripts/15_main_figures.R
# ======================================================================

source("scripts/utils_functions.R")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(patchwork); library(png); library(grid); library(ragg)
})

TAB <- P_TAB
OUT <- P_FIG

COL_COMP <- "#2a78d6"   # composition (blue)
COL_ABUN <- "#eb6834"   # abundance (orange)
COL_REF  <- "#8a8984"   # baselines / reference marks
INK      <- "#0b0b0b"
INK2     <- "#52514e"

theme_fig <- function() {
  theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = INK),
          axis.title = element_text(colour = INK),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          panel.grid.major.x = element_line(colour = "#ecebe7", linewidth = 0.3),
          plot.tag = element_text(face = "bold", size = 16),
          legend.position = "none")
}

taxa <- c(D11 = "Pedobacter sp.", D14 = "Microbacterium sp.",
          D17 = "Curtobacterium sp.", I2 = "Yersinia sp.",
          I8 = "P. fluorescens", I9 = "Chromobacterium sp.",
          I11 = "P. anguilliseptica", I15 = "Serratia sp.",
          I18 = "Aeromonas sp.", I20 = "Chryseobacterium sp.",
          I22 = "Erwinia sp.", I23 = "Janthinobacterium sp.")

# ---------------------------------------------------------------------
# Figure 4: existing Fig 3 (a, b) + decomposition ladder (c) + WTA (d)
# ---------------------------------------------------------------------
s11 <- read_csv(TAB("Table_S11_model_overall_metrics.csv"), show_col_types = FALSE)
s19 <- read_csv(TAB("Table_S19_biomass_model_decomposition.csv"), show_col_types = FALSE)
s25 <- read_csv(TAB("Table_S25_winner_takes_all.csv"), show_col_types = FALSE)
s22 <- read_csv(TAB("Table_S22_trait_comparison.csv"), show_col_types = FALSE)

ladder <- tibble::tibble(
  model = c("Stress + richness only (no trait)",
            "Growth rates, equal shares",
            "Growth rates, winner-takes-all",
            "Growth rates, observed composition (oracle)",
            "Growth rates, predicted composition (model)"),
  R2 = c(s19$R2[grepl("^Null", s19$Model)],
         s11$R2[grepl("^Equal", s11$Model)],
         s25$OD_R2,
         s11$R2[grepl("oracle", s11$Model)],
         s19$R2[grepl("^Full", s19$Model)]),
  is_model = c(FALSE, FALSE, FALSE, FALSE, TRUE)
) %>% mutate(model = factor(model, levels = rev(model)))

p4c <- ggplot(ladder, aes(x = R2, y = model)) +
  geom_segment(aes(x = 0, xend = R2, yend = model), colour = COL_REF, linewidth = 0.6) +
  geom_point(aes(colour = is_model), size = 3.2) +
  geom_text(aes(label = sprintf("%.3f", R2)), hjust = -0.35, size = 3.4, colour = INK2) +
  scale_colour_manual(values = c(`FALSE` = COL_REF, `TRUE` = COL_ABUN)) +
  scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = expression("Community abundance "*R^2*" (OD"[600]*")"), y = NULL, tag = "c") +
  theme_fig()

wta <- tibble::tibble(
  response = factor(c("Composition\n(weighted R²)", "Composition\n(weighted R²)",
                      "Abundance\n(R²)", "Abundance\n(R²)"),
                    levels = c("Abundance\n(R²)", "Composition\n(weighted R²)")),
  predictor = c("Winner-takes-all", "Model", "Winner-takes-all", "Model"),
  value = c(s25$comp_wR2, s22$comp_wR2[s22$trait == "r_gompertz"], s25$OD_R2, s19$R2[grepl("^Full", s19$Model)])
)
wta <- wta %>% mutate(y = as.numeric(response) + ifelse(predictor == "Model", 0.13, -0.13))
p4d <- ggplot(wta, aes(x = value, y = y)) +
  geom_vline(xintercept = 0, colour = COL_REF, linewidth = 0.4, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = value, yend = y), colour = COL_REF, linewidth = 0.5) +
  geom_point(aes(shape = predictor, fill = predictor), size = 3.4, colour = "white", stroke = 0.8) +
  geom_text(aes(label = sprintf("%.3f", value), hjust = ifelse(value < 0, 1.35, -0.35)),
            size = 3.2, colour = INK2) +
  scale_shape_manual(values = c(Model = 21, `Winner-takes-all` = 24)) +
  scale_fill_manual(values = c(Model = COL_COMP, `Winner-takes-all` = COL_ABUN)) +
  scale_x_continuous(limits = c(-1.3, 1.2), breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(breaks = seq_along(levels(wta$response)), labels = levels(wta$response),
                     limits = c(0.5, 2.5)) +
  labs(x = expression(R^2), y = NULL, tag = "d", shape = NULL, fill = NULL) +
  theme_fig() + theme(legend.position = "top", legend.justification = "left",
                      legend.margin = margin(0, 0, 0, 0))

top <- wrap_elements(full = rasterGrob(readPNG(P_FIG("Fig_3.png")), interpolate = TRUE))
fig4 <- top / (p4c | p4d + plot_layout(widths = 1)) + plot_layout(heights = c(1.05, 0.75))
ggsave(OUT("Fig_4.png"), fig4, width = 260, height = 175, units = "mm", dpi = 300, device = ragg::agg_png, bg = "white")

# ---------------------------------------------------------------------
# Figure 5: transfer
# ---------------------------------------------------------------------
s15 <- read_csv(TAB("Table_S15_blockedcv_com_id_bias1_by_fold.csv"), show_col_types = FALSE)
s16 <- read_csv(TAB("Table_S16_blockedcv_com_id_bias1_by_stress.csv"), show_col_types = FALSE)
s20 <- read_csv(TAB("Table_S20_losocv_by_stress.csv"), show_col_types = FALSE)
s28 <- read_csv(TAB("Table_S28_singles_to_combos.csv"), show_col_types = FALSE)
s27 <- read_csv(TAB("Table_S27_leave_one_taxon_out.csv"), show_col_types = FALSE)

tests <- c("Unseen\ncommunities", "Unseen stress\nregime", "Combinations\nfrom singles")

comp <- bind_rows(
  tibble::tibble(test = tests[1], value = s16$RMSE),
  tibble::tibble(test = tests[2], value = s20$comp_RMSE),
  tibble::tibble(test = tests[3], value = s28$comp_RMSE[!grepl("POOLED", s28$Stress)])
) %>% mutate(test = factor(test, levels = tests))

abun <- bind_rows(
  tibble::tibble(test = tests[1], value = s15$od_R2, pooled = FALSE),
  tibble::tibble(test = tests[2], value = s20$OD_R2, pooled = FALSE),
  tibble::tibble(test = tests[3], value = s28$OD_R2_pooled[grepl("POOLED", s28$Stress)], pooled = TRUE)
) %>% mutate(test = factor(test, levels = tests))

set.seed(1)
p5a <- ggplot(comp, aes(x = test, y = value)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = INK2) +
  geom_point(position = position_jitter(width = 0.12, height = 0), colour = COL_COMP,
             size = 2.6, alpha = 0.9) +
  scale_y_continuous(limits = c(0, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Held-out composition RMSE", tag = "a",
       title = "Composition transfers") +
  theme_fig() + theme(panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(colour = "#ecebe7", linewidth = 0.3),
                      plot.title = element_text(size = 11, face = "bold"))

p5b <- ggplot(abun, aes(x = test, y = value)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_hline(yintercept = s19$R2[grepl("^Full", s19$Model)], colour = COL_REF,
             linewidth = 0.4, linetype = "dotted") +
  geom_point(data = filter(abun, !pooled), position = position_jitter(width = 0.12, height = 0),
             colour = COL_ABUN, size = 2.6, alpha = 0.9) +
  geom_point(data = filter(abun, pooled), shape = 23, fill = COL_ABUN, colour = "white",
             size = 3.6, stroke = 0.8) +
  annotate("text", x = 3.45, y = 0.88, label = "in-sample", hjust = 1, vjust = -0.5,
           size = 3, colour = INK2) +
  scale_y_continuous(limits = c(-11, 1.2), breaks = seq(-10, 0, 2.5)) +
  labs(x = NULL, y = expression("Held-out abundance "*R^2), tag = "b",
       title = "Absolute abundance does not") +
  theme_fig() + theme(panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(colour = "#ecebe7", linewidth = 0.3),
                      plot.title = element_text(size = 11, face = "bold"))

loto <- s27 %>% filter(!grepl("POOLED", taxon)) %>%
  mutate(name = taxa[taxon], name = factor(name, levels = name[order(focal_spearman)]))
pooled <- s27 %>% filter(grepl("POOLED", taxon))

p5c <- ggplot(loto, aes(x = focal_spearman, y = name)) +
  geom_vline(xintercept = pooled$focal_spearman, colour = COL_REF, linetype = "dotted", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = focal_spearman, yend = name), colour = COL_REF, linewidth = 0.5) +
  geom_point(colour = COL_COMP, size = 2.8) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = expression("Spearman "*rho*", held-out taxon"), y = NULL, tag = "c") +
  theme_fig() + theme(axis.text.y = element_text(face = "italic"))

loto_long <- bind_rows(
  loto %>% transmute(name, role = "Held-out taxon", value = focal_RMSE),
  loto %>% transmute(name, role = "Calibrated taxa", value = others_RMSE)
)
p5d <- ggplot(loto_long, aes(x = value, y = name)) +
  geom_line(aes(group = name), colour = COL_REF, linewidth = 0.5) +
  geom_point(aes(fill = role, shape = role), size = 2.8, colour = "white", stroke = 0.6) +
  scale_fill_manual(values = c(`Held-out taxon` = COL_COMP, `Calibrated taxa` = COL_REF)) +
  scale_shape_manual(values = c(`Held-out taxon` = 21, `Calibrated taxa` = 22)) +
  scale_x_continuous(limits = c(0, 0.35), expand = c(0, 0)) +
  labs(x = "Relative-abundance RMSE", y = NULL, tag = "d", fill = NULL, shape = NULL) +
  theme_fig() + theme(axis.text.y = element_blank(), legend.position = "top",
                      legend.justification = "left", legend.margin = margin(0, 0, 0, 0))

fig5 <- (p5a | p5b) / (p5c | p5d) + plot_layout(heights = c(1, 1.1))
ggsave(OUT("Fig_5.png"), fig5, width = 200, height = 190, units = "mm", dpi = 300, device = ragg::agg_png, bg = "white")

message("Wrote ", P_FIG("Fig_4.png"), " and ", P_FIG("Fig_5.png"))
