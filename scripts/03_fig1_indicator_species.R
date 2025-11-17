# ======================================================================
# Fig 1: Indicator Species Analysis Across Stress Treatments
#
# Products:
# - Fig 1d: Indicator taxa (Top 1 per stress condition) plotted by Genus
# - Table S6: Full indicator species statistics for all ASVs × stress groups
#
# Overview:
# - Loads cleaned ASV abundance matrix, aligned metadata, and taxonomy table
# - Identifies which ASVs are most indicative of each stress condition
#   using indicator value analysis (IndVal.g)
# - Cleans and harmonises genus-level taxonomy (including Pseudomonas variants)
# - Extracts the top-N indicator taxa per stress condition
# - Saves the full indicator statistics table and generates Fig 1d
#
# Folder conventions:
# - Inputs (cleaned):  rds/       (via P_RDS())
# - Plots:             plots/     (via P_FIG())
# - Tables (CSV):      tables/    (via P_TAB())
#
# Outline:
# 1) Load ASV matrix, metadata, and taxonomy
# 2) Determine which stress conditions are present and map to pretty labels
# 3) Run indicator species analysis (IndVal.g) for each stress vs. rest
# 4) Combine all indicator stats, select Top N per condition
# 5) Save full indicator stats as Table S6
# 6) Plot Top N indicator genera by stress condition and save Fig 1d

# ======================================================================
#  Fig 1 (indicator species) -> Table S6, Fig 1d
source("scripts/utils_functions.R"); ensure_packages()

ASV_t_clean        <- readRDS(P_RDS("ASV_matrix_clean.rds"))
env_option_A_clean <- readRDS(P_RDS("metadata_clean.rds"))
taxa_table         <- readRDS(P_RDS("taxa_table.rds"))

grp <- env_option_A_clean$Stress
stopifnot(nrow(ASV_t_clean) == length(grp))

short_levels   <- names(pretty_map)
present_short  <- short_levels[short_levels %in% unique(as.character(grp))]
present_pretty <- unname(pretty_map[present_short])

tax_df <- as.data.frame(taxa_table) %>% tibble::rownames_to_column("ASV")

TopN       <- 1
assoc_func <- "IndVal.g"

panel_list <- lapply(present_short, function(g) {
  grp_bin <- factor(ifelse(as.character(grp) == g, g, "Other"), levels = c("Other", g))
  fit <- indicspecies::multipatt(ASV_t_clean, grp_bin, func = assoc_func, max.order = 1, duleg = FALSE,
                                 control = permute::how(nperm = 0))
  sig <- as.data.frame(fit$sign); sig$ASV <- rownames(sig)
  mem_col <- paste0("s.", levels(grp_bin)[2])
  sig %>% dplyr::transmute(ASV, stat = as.numeric(stat), mem = .data[[mem_col]]) %>%
    dplyr::filter(is.finite(stat), stat > 0, mem == 1) %>% dplyr::select(-mem) %>%
    dplyr::left_join(tax_df, by = "ASV") %>%
    dplyr::mutate(
      Genus = stringr::str_remove(coalesce(Genus, ""), "^g__"),
      Genus = ifelse(Genus == "" | is.na(Genus), ASV, Genus),
      Genus = dplyr::case_when(Genus == "Pseudomonas"  ~ "Pseudomonas.A",
                               Genus == "Pseudomonas2" ~ "Pseudomonas.F",
                               TRUE ~ Genus),
      Group_short  = g,
      Group_pretty = pretty_map[g]
    ) %>% dplyr::arrange(dplyr::desc(stat))
})

all_isa_stats <- dplyr::bind_rows(panel_list) %>%
  dplyr::mutate(Group_pretty = factor(Group_pretty, levels = present_pretty))

topN_by_group <- all_isa_stats %>% dplyr::group_by(Group_short) %>% dplyr::slice_head(n = TopN) %>% dplyr::ungroup()

# Save Table S6 (rounded to 3 decimals) and drop Group_pretty
all_isa_stats_out <- all_isa_stats %>% dplyr::select(-Group_pretty) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))
readr::write_csv(all_isa_stats_out, P_TAB("Table_S6_ISA_stats_full.csv"))

# Plot
rng   <- range(topN_by_group$stat, na.rm = TRUE)
brks0 <- pretty(rng, n = 3); brks  <- if (length(brks0) >= 3) brks0[c(1, floor(length(brks0)/2), length(brks0))] else rng
labf  <- scales::number_format(accuracy = 0.01)
col_vec <- c("#e8f1fb", "#6aa0d8", "#08306b")
POINT_SIZE <- 5.5

p_isa <- ggplot(topN_by_group, aes(x = Group_pretty, y = Genus)) +
  geom_point(aes(colour = stat), size = POINT_SIZE, shape = 16, alpha = 0.9) +
  scale_colour_gradientn(colours = col_vec, limits  = rng, breaks  = brks, labels  = labf, name = "Indicator value") +
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5, label.position = "bottom",
                                  barwidth = unit(80, "pt"), barheight = unit(6, "pt"), ticks.colour = "grey20")) +
  scale_x_discrete(labels = scales::label_wrap(14)) +
  labs(x = "Stress condition", y = "Indicator taxon (Genus)") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(legend.position   = "top", legend.background = element_rect(fill = "white", colour = NA),
        legend.key        = element_rect(fill = "white", colour = NA), legend.title = element_text(size = 14, face = "bold"),
        legend.text       = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y       = element_text(size = 14, face = "italic"), axis.title.x = element_text(size = 16),
        axis.title.y      = element_text(size = 16), axis.line = element_line(linewidth = 0.4),
        axis.ticks.length = unit(3, "pt"), plot.margin = margin(8, 12, 8, 8))

ggsave(P_FIG("Fig_1d_ISA_Top1_1vsRest_singlepanel.tiff"), label_plot(p_isa, "1d"),
       width = 18, height = 15, units = "cm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")
