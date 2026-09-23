# ======================================================================
# 16_fig1_schematic.R
#
# Figure 1: design schematic. No data. Drawn to the specification in the
# comment above @fig-1 in manuscript/manuscript.qmd:
#
#   a  the 2 x 2 x 2 factorial, the taxon pool and the 22 communities
#   b  step 1: monoculture reaction norms give a taxon x regime table of
#      growth rates; rank order changes between control and stress level
#   c  step 2: growth rates -> softmax -> composition -> AWM -> abundance
#   d  transfer map: what was withheld, and what carried over
#
# Colours follow the rest of the paper: stress regimes from pal_short
# (utils_functions.R), composition in COL_COMP and abundance in COL_ABUN
# (15_main_figures.R).
#
# Run from the repository root:  Rscript scripts/16_fig1_schematic.R
# Writes results/figures/Fig_1_design.png and .pdf
# ======================================================================

source("scripts/utils_functions.R")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(grid); library(ragg)
})

COL_COMP <- "#2a78d6"; COL_ABUN <- "#eb6834"; COL_REF <- "#8a8984"
INK <- "#0b0b0b"; INK2 <- "#52514e"; PAPER <- "#f6f5f2"
BASE <- 9
# Symbols (arrows, degree, multiplication) need a UTF-8 session locale; on
# a bare Linux box run with  LC_ALL=C.UTF-8 Rscript scripts/16_fig1_schematic.R

theme_blank <- function() {
  theme_void(base_size = BASE) +
    theme(plot.tag = element_text(face = "bold", size = 16),
          plot.title = element_text(size = BASE + 1, face = "bold", colour = INK, hjust = 0),
          plot.margin = margin(4, 4, 4, 4))
}
box <- function(x0, y0, x1, y1, fill = "white", col = COL_REF, lwd = 0.4, r = 0.02) {
  annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = col, linewidth = lwd)
}
arr <- function(x0, y0, x1, y1, col = INK2, lwd = 0.5, len = 5) {
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1, colour = col, linewidth = lwd,
           arrow = arrow(length = unit(len, "pt"), type = "closed"), linejoin = "mitre")
}
txt <- function(x, y, label, size = 2.9, col = INK, ...) annotate("text", x = x, y = y, label = label, size = size, colour = col, ...)
tick <- function(x, y, col = COL_COMP, s = 0.045) {
  list(annotate("segment", x = x - s, y = y, xend = x - s * 0.25, yend = y - s * 0.8, colour = col, linewidth = 0.9, lineend = "round"),
       annotate("segment", x = x - s * 0.25, y = y - s * 0.8, xend = x + s * 1.1, yend = y + s * 0.9, colour = col, linewidth = 0.9, lineend = "round"))
}
cross <- function(x, y, col = COL_ABUN, s = 0.04) {
  list(annotate("segment", x = x - s, y = y - s, xend = x + s, yend = y + s, colour = col, linewidth = 0.9, lineend = "round"),
       annotate("segment", x = x - s, y = y + s, xend = x + s, yend = y - s, colour = col, linewidth = 0.9, lineend = "round"))
}

# ======================================================================
# a: design
# ======================================================================
# isometric cube: x = temperature, y = salinity, z = pH (into the page)
proj <- function(x, y, z) c(X = x + 0.45 * z, Y = y + 0.32 * z)
corners <- tribble(
  ~T, ~S, ~P, ~regime,
  0, 0, 0, "Control",  1, 0, 0, "Temp",     0, 1, 0, "Sal",      0, 0, 1, "pH",
  1, 1, 0, "SalTemp",  1, 0, 1, "pHTemp",   0, 1, 1, "pHSal",    1, 1, 1, "pHSalTemp"
) %>% rowwise() %>% mutate(X = proj(T, S, P)["X"], Y = proj(T, S, P)["Y"]) %>% ungroup() %>%
  mutate(lab = pretty_map[regime], lab = gsub(" x ", " × ", lab),
         lab = case_when(regime == "pHSalTemp" ~ "pH × Sal × Temp",
                         regime == "SalTemp" ~ "Sal × Temp", regime == "pHTemp" ~ "pH × Temp",
                         regime == "pHSal" ~ "pH × Sal", TRUE ~ regime),
         hj = ifelse(T == 1, -0.18, 1.18), vj = ifelse(S == 1, -0.55, 1.55))
edges <- bind_rows(lapply(seq_len(nrow(corners)), function(i) {
  a <- corners[i, ]
  bind_rows(lapply(seq_len(nrow(corners)), function(j) {
    b <- corners[j, ]
    if (j > i && sum(abs(c(a$T - b$T, a$S - b$S, a$P - b$P))) == 1)
      tibble(x = a$X, y = a$Y, xend = b$X, yend = b$Y, back = (a$P + b$P) == 2 || (a$T == 0 & b$T == 0 & a$P + b$P >= 1 & a$S + b$S == 2 * a$S))
  }))
}))
cube_off <- c(0.42, 0.42); cube_sc <- 0.52
cx <- function(v) cube_off[1] + v * cube_sc; cy <- function(v) cube_off[2] + v * cube_sc

pa <- ggplot() +
  annotate("segment", x = cx(edges$x), y = cy(edges$y), xend = cx(edges$xend), yend = cy(edges$yend),
           colour = COL_REF, linewidth = 0.45) +
  geom_point(data = corners, aes(cx(X), cy(Y), fill = regime), shape = 21, size = 4.6, colour = "white", stroke = 0.7) +
  geom_text(data = corners, aes(cx(X), cy(Y), label = lab, hjust = hj, vjust = vj), size = 2.5, colour = INK) +
  scale_fill_manual(values = pal_short, guide = "none") +
  # axis arrows
  arr(cx(-0.02), cy(-0.30), cx(1.02), cy(-0.30), col = INK2, lwd = 0.45, len = 4) +
  txt(cx(0.5), cy(-0.40), "Temperature: 20 → 38 °C", size = 2.5, col = INK2) +
  arr(cx(-0.56), cy(-0.02), cx(-0.56), cy(1.02), col = INK2, lwd = 0.45, len = 4) +
  txt(cx(-0.68), cy(0.5), "Salinity: 0 → 20 g NaCl/L", size = 2.5, col = INK2, angle = 90) +
  arr(cx(1.30), cy(-0.30), cx(1.30 + 0.45), cy(-0.30 + 0.32), col = INK2, lwd = 0.45, len = 4) +
  txt(cx(1.82), cy(-0.02), "pH: 7.2 → 5.5", size = 2.5, col = INK2, hjust = 0) +
  txt(0.05, 1.36, "Eight stress regimes, 2 × 2 × 2", size = 3, col = INK, fontface = "bold", hjust = 0) +
  # taxon pool
  txt(2.55, 1.36, "Twelve taxa", size = 3, col = INK, fontface = "bold", hjust = 0) +
  geom_point(data = expand.grid(i = 0:3, j = 0:2) %>% mutate(x = 2.6 + i * 0.13, y = 1.2 - j * 0.13),
             aes(x, y), shape = 21, size = 3.2, fill = "#9aa5b1", colour = "white", stroke = 0.6) +
  arr(2.83, 0.86, 2.83, 0.72, len = 4) +
  txt(2.83, 0.66, "assembled into 22 compositions", size = 2.5, col = INK2) +
  # community boxes
  box(2.20, 0.18, 2.58, 0.52) + box(2.64, 0.18, 3.02, 0.52) + box(3.08, 0.18, 3.46, 0.52) +
  geom_point(data = tibble(x = c(2.32, 2.46), y = 0.40), aes(x, y), shape = 21, size = 2.6, fill = "#9aa5b1", colour = "white", stroke = 0.5) +
  geom_point(data = expand.grid(i = 0:1, j = 0:1) %>% mutate(x = 2.76 + i * 0.14, y = 0.45 - j * 0.11),
             aes(x, y), shape = 21, size = 2.6, fill = "#9aa5b1", colour = "white", stroke = 0.5) +
  geom_point(data = expand.grid(i = 0:3, j = 0:1) %>% mutate(x = 3.15 + i * 0.08, y = 0.45 - j * 0.11),
             aes(x, y), shape = 21, size = 2.6, fill = "#9aa5b1", colour = "white", stroke = 0.5) +
  txt(2.39, 0.25, "2 taxa (6)", size = 2.3, col = INK2) + txt(2.83, 0.25, "4 taxa (6)", size = 2.3, col = INK2) +
  txt(3.27, 0.25, "8 taxa (10)", size = 2.3, col = INK2) +
  # every community in every regime
  arr(2.12, 0.24, 1.78, 0.24, len = 4) +
  txt(1.95, 0.13, "each grown in\nevery regime", size = 2.4, col = INK2, lineheight = 0.9) +
  txt(2.83, -0.04, "48 h batch culture → 16S amplicon composition + OD₆₀₀", size = 2.5, col = INK) +
  coord_equal(xlim = c(-0.02, 3.55), ylim = c(-0.10, 1.42), expand = FALSE, clip = "off") +
  theme_blank()

# ======================================================================
# b: step 1, reaction norms
# ======================================================================
xs <- seq(0, 1, length.out = 200)
norms <- bind_rows(
  tibble(taxon = "A", x = xs, y = 1.00 * exp(-((xs - 0.20) / 0.28)^2)),
  tibble(taxon = "B", x = xs, y = 0.85 * exp(-((xs - 0.45) / 0.30)^2)),
  tibble(taxon = "C", x = xs, y = 0.70 * exp(-((xs - 0.70) / 0.28)^2)),
  tibble(taxon = "D", x = xs, y = 0.55 * exp(-((xs - 0.85) / 0.30)^2)))
tax_cols <- c(A = "#4f5d75", B = "#8d99ae", C = "#bc6c25", D = "#606c38")
x_ctrl <- 0.22; x_str <- 0.78
at <- norms %>% filter(abs(x - x_ctrl) < 0.003 | abs(x - x_str) < 0.003) %>%
  mutate(level = ifelse(x < 0.5, "control", "stress")) %>% group_by(level) %>%
  mutate(rank = rank(-y)) %>% ungroup()

pb <- ggplot() +
  annotate("segment", x = c(x_ctrl, x_str), xend = c(x_ctrl, x_str), y = 0, yend = 1.08,
           colour = COL_REF, linetype = "22", linewidth = 0.4) +
  geom_line(data = norms, aes(x, y, colour = taxon), linewidth = 0.9) +
  geom_point(data = at, aes(x, y, fill = taxon), shape = 21, size = 2.8, colour = "white", stroke = 0.6) +
  scale_colour_manual(values = tax_cols, guide = "none") + scale_fill_manual(values = tax_cols, guide = "none") +
  txt(x_ctrl, 1.15, "control\nlevel", size = 2.5, col = INK2, lineheight = 0.9) +
  txt(x_str, 1.15, "stress\nlevel", size = 2.5, col = INK2, lineheight = 0.9) +
  # rank labels at each level
  geom_text(data = at %>% filter(level == "control") %>% mutate(y = ifelse(taxon == "D", y - 0.055, ifelse(taxon == "C", y + 0.05, y))),
            aes(x = x - 0.04, y = y, label = taxon), size = 2.4, colour = INK2, hjust = 1) +
  geom_text(data = at %>% filter(level == "stress"), aes(x = x + 0.04, y = y, label = taxon), size = 2.4, colour = INK2, hjust = 0) +
  txt(0.5, -0.16, "stressor level", size = 2.8, col = INK) +
  txt(-0.12, 0.5, "monoculture growth rate", size = 2.8, col = INK, angle = 90) +
  annotate("segment", x = 0, xend = 1.02, y = 0, yend = 0, colour = INK2, linewidth = 0.4) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 1.08, colour = INK2, linewidth = 0.4) +
  txt(0.5, 1.42, "Reaction norms measured once per taxon (Carmichael et al. 2025)", size = 2.7, col = INK2) +
  txt(0.5, 1.32, "rank order changes between levels", size = 2.6, col = INK, fontface = "italic") +
  # output: taxon x regime table
  box(1.22, 0.15, 1.95, 0.95, fill = PAPER) +
  txt(1.585, 0.86, "growth rate  g", size = 2.7, col = INK, fontface = "bold") +
  txt(1.585, 0.76, "taxon × regime", size = 2.4, col = INK2) +
  geom_tile(data = { set.seed(1); expand.grid(i = 1:8, j = 1:6) %>% mutate(x = 1.29 + (i - 1) * 0.082, y = 0.63 - (j - 1) * 0.078,
                                                            v = runif(48, 0.15, 1)) },
            aes(x, y, fill = NULL, alpha = v), fill = COL_COMP, width = 0.072, height = 0.068) +
  scale_alpha_identity() +
  txt(1.585, 0.19, "8 regimes × 12 taxa", size = 2.2, col = INK2) +
  arr(1.06, 0.55, 1.19, 0.55, len = 4) +
  coord_cartesian(xlim = c(-0.18, 1.98), ylim = c(-0.22, 1.5), expand = FALSE, clip = "off") +
  theme_blank()

# ======================================================================
# c: step 2, prediction
# ======================================================================
stack <- tibble(taxon = c("A", "B", "C", "D"), p = c(0.46, 0.28, 0.17, 0.09)) %>%
  mutate(top = cumsum(p), bottom = top - p)
blues <- c(A = "#1f5fae", B = "#2a78d6", C = "#6ea3e6", D = "#b3cff2")
pts <- tibble(x = c(0.15, 0.32, 0.48, 0.62, 0.80, 0.90), y = c(0.22, 0.36, 0.42, 0.58, 0.70, 0.84))

pc <- ggplot() +
  # 1 growth rates
  box(0.00, 0.30, 0.62, 0.90, fill = PAPER) +
  txt(0.31, 0.79, "growth rates", size = 2.7, col = INK, fontface = "bold") +
  txt(0.31, 0.65, "g, z-scored\nwithin regime", size = 2.4, col = INK2, lineheight = 0.9) +
  txt(0.31, 0.43, "+ taxon offsets δ", size = 2.4, col = INK2) +
  arr(0.65, 0.60, 0.93, 0.60, len = 4) + txt(0.79, 0.97, "softmax", size = 2.3, col = INK2) +
  # 2 composition
  box(0.96, 0.30, 1.58, 0.90, fill = "white") +
  txt(1.27, 0.79, "composition", size = 2.7, col = COL_COMP, fontface = "bold") +
  geom_rect(data = stack, aes(xmin = 1.13, xmax = 1.41, ymin = 0.37 + bottom * 0.33, ymax = 0.37 + top * 0.33, fill = taxon), colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = blues, guide = "none") +
  txt(1.27, 0.335, "predicted relative abundance", size = 1.9, col = INK2) +
  arr(1.61, 0.60, 1.89, 0.60, len = 4) + txt(1.75, 0.97, "weight by g", size = 2.3, col = INK2) +
  # 3 AWM
  box(1.92, 0.30, 2.54, 0.90, fill = PAPER) +
  txt(2.23, 0.79, "AWM", size = 2.7, col = INK, fontface = "bold") +
  txt(2.23, 0.63, "abundance-weighted\nmean growth", size = 2.3, col = INK2, lineheight = 0.9) +
  txt(2.23, 0.43, "Σ p · g", size = 2.8, col = INK2) +
  arr(2.57, 0.60, 2.85, 0.60, len = 4) + txt(2.71, 0.97, "regression", size = 2.3, col = INK2) +
  txt(2.71, 0.47, "regime-\nspecific", size = 2.1, col = COL_ABUN, lineheight = 0.9) +
  # 4 abundance
  box(2.88, 0.30, 3.50, 0.90, fill = "white") +
  txt(3.19, 0.79, "abundance", size = 2.7, col = COL_ABUN, fontface = "bold") +
  geom_point(data = pts, aes(2.97 + x * 0.44, 0.37 + y * 0.33), colour = COL_ABUN, size = 1.4) +
  annotate("segment", x = 3.00, xend = 3.38, y = 0.37 + 0.14 * 0.33, yend = 0.37 + 0.90 * 0.33, colour = COL_ABUN, linewidth = 0.6) +
  txt(3.19, 0.335, "predicted OD₆₀₀", size = 1.9, col = INK2) +
  # brace labels
  annotate("segment", x = 0.00, xend = 1.58, y = 0.17, yend = 0.17, colour = COL_COMP, linewidth = 0.8) +
  txt(0.79, 0.07, "transfers without community data from the focal regime", size = 2.3, col = COL_COMP) +
  annotate("segment", x = 1.92, xend = 3.50, y = 0.17, yend = 0.17, colour = COL_ABUN, linewidth = 0.8) +
  txt(2.71, 0.07, "needs calibration against the focal regime", size = 2.3, col = COL_ABUN) +
  coord_equal(xlim = c(-0.05, 3.55), ylim = c(0.0, 1.04), expand = FALSE, clip = "off") +
  theme_blank()

# ======================================================================
# d: transfer map
# ======================================================================
tiles <- tibble(
  i = 1:4,
  title = c("Unseen\ncommunities", "Unseen\ntaxa", "Unseen\nstress regimes", "Combinations\nfrom single stressors"),
  held = c("22 compositions,\nblocked 5-fold", "each taxon absent\nfrom every fitted\ncommunity", "each regime absent\nfrom fitting", "fitted on single\nstressors only"),
  comp = c(TRUE, TRUE, TRUE, TRUE),
  abun = c("yes", "na", "no", "no")) %>%
  mutate(x0 = (i - 1) * 0.9, x1 = x0 + 0.82)

pd <- ggplot() + theme_blank()
for (k in seq_len(nrow(tiles))) {
  t <- tiles[k, ]; xm <- (t$x0 + t$x1) / 2
  pd <- pd + box(t$x0, 0.0, t$x1, 1.0, fill = PAPER, col = COL_REF) +
    txt(xm, 0.88, t$title, size = 2.7, col = INK, fontface = "bold", lineheight = 0.9) +
    txt(xm, 0.66, t$held, size = 2.2, col = INK2, lineheight = 0.9) +
    annotate("segment", x = t$x0 + 0.06, xend = t$x1 - 0.06, y = 0.47, yend = 0.47, colour = COL_REF, linewidth = 0.3) +
    txt(t$x0 + 0.08, 0.36, "composition", size = 2.4, col = COL_COMP, hjust = 0) +
    txt(t$x0 + 0.08, 0.16, "abundance", size = 2.4, col = COL_ABUN, hjust = 0) +
    tick(t$x1 - 0.13, 0.36, col = COL_COMP)
  pd <- pd + switch(t$abun,
    yes = tick(t$x1 - 0.13, 0.16, col = COL_ABUN),
    no  = cross(t$x1 - 0.13, 0.16, col = COL_ABUN),
    na  = list(txt(t$x1 - 0.13, 0.16, "–", size = 3, col = COL_REF)))
}
pd <- pd +
  txt(1.76, 1.16, "Held-out test", size = 2.8, col = INK, fontface = "bold") +
  tick(0.62, -0.12, col = COL_COMP, s = 0.035) + txt(0.70, -0.12, "transfers", size = 2.3, col = INK2, hjust = 0) +
  cross(1.42, -0.12, col = COL_ABUN, s = 0.032) + txt(1.50, -0.12, "requires focal-regime calibration", size = 2.3, col = INK2, hjust = 0) +
  txt(2.92, -0.12, "–", size = 3, col = COL_REF) + txt(3.00, -0.12, "not assessed", size = 2.3, col = INK2, hjust = 0) +
  coord_equal(xlim = c(-0.05, 3.57), ylim = c(-0.2, 1.25), expand = FALSE, clip = "off")

# ======================================================================
# assemble
# ======================================================================
fig <- (pa + labs(tag = "a")) / (pb + labs(tag = "b")) / (pc + labs(tag = "c")) / (pd + labs(tag = "d")) +
  plot_layout(heights = c(1.52, 1.6, 1.04, 1.45))

ggsave(P_FIG("Fig_1_design.png"), fig, width = 180, height = 235, units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(P_FIG("Fig_1_design.pdf"), fig, width = 180, height = 235, units = "mm", bg = "white", device = cairo_pdf)
message("Done: results/figures/Fig_1_design.png and .pdf")
