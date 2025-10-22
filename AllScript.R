# Commmunity analysis - linear model - 4/12/23

#Step 1: Clear the Brain (consloe)
rm(list = ls())

#load libraries
library("ggplot2")
library("gridExtra")
require(reshape2)
library("readxl")
library(tidyverse)
library(ggpubr)
library(broom)
library(rstatix)
library(car)
library(emmeans)
library(AICcmodavg)
library(lmtest)

# Get the data
All_data <- read.csv('CommunityOD_All_Jan23_edit_zeros.csv')

# Mean center the data - create function
center_scale <- function(x) {
  scale(x, center = TRUE, scale = FALSE)
}

#log OD and only look at TF and mean centre the data
All_data <- All_data %>% mutate(log10_OD = log10(OD))%>% filter(t == 2) 

#remove evo data 
All_data_noEvo <- All_data %>% filter(Evo_Treatment != "Evo") %>%
  select(Id, Stress, TempLevel, Community, Media, Diversity_Level, OD, log10_OD)

# calculate the mean functioing per monoculture
All_data_meanMono <- All_data_noEvo %>% filter(Diversity_Level == 1) %>%
  group_by(Id, Stress, TempLevel, Community, Media, Diversity_Level) %>%
  summarise(OD = mean(OD), 
            log10_OD = mean(log10_OD))

All_data_noEvo <- All_data_noEvo %>% filter(Diversity_Level != 1)

# merge the two
All_data_noEvo <- rbind(All_data_noEvo, All_data_meanMono)

#calculate average biomass for each diversity
mean_Biomass <- All_data_noEvo %>% group_by(Diversity_Level, Stress) %>%
  summarise(log_mean_biomass = mean(log10_OD)) %>% #log first then average
  filter(complete.cases(Diversity_Level)) %>%
  mutate(log10_Diversity_centered = log10(Diversity_Level)) 

# calculate centered diversity
median <- median(mean_Biomass$log10_Diversity_centered)

mean_Biomass <- mean_Biomass %>%
  mutate(log10_Diversity_centered = log10(Diversity_Level) - median) 

All_data_noEvo <- All_data_noEvo %>%
  mutate(log10_Diversity_centered = log10(Diversity_Level) - median) 


write.csv(All_data_noEvo, "Community_BEF_data.csv", row.names = FALSE)


# view general patterns
quickPlot <- ggplot(mean_Biomass, aes(log10_Diversity_centered,log_mean_biomass )) +
  geom_point()+
  facet_grid(~Stress) +
  facet_wrap(~fct_relevel(Stress, "Control", "pH", "Sal", "Temp", "pH_Sal", "pH_Temp", "Sal_Temp", "pH_Sal_Temp"), nrow = 1)
quickPlot

#### linear model with all data ####
Diversity_mod1 <- lm(log10_OD ~ as.vector(log10_Diversity_centered) * Stress, data = All_data_noEvo )
Diversity_mod2 <- lm(log10_OD ~ as.vector(log10_Diversity_centered) + Stress, data = All_data_noEvo )
Diversity_mod3 <- lm(log10_OD ~ as.vector(log10_Diversity_centered), data = All_data_noEvo )
Diversity_mod4 <- lm(log10_OD ~Stress, data = All_data_noEvo )

Three_way <- All_data_noEvo  %>% filter(Stress == "pH_Sal_Temp")
#Diversity_mod1 <- lm(log10_OD ~ as.vector(log10_Diversity_centered) , data = Three_way )

#plot(Diversity_mod1)
AICc(Diversity_mod1)
AICc(Diversity_mod2)
AICc(Diversity_mod3)
AICc(Diversity_mod4)

lrtest(Diversity_mod1, Diversity_mod2) 

anova(Diversity_mod1) # significant effects of diversity level on functioning. The Stress treatment also significantly
# effects the functioning (changes in the intercept) and the way in which diversity level effects functioning varies by stress treatment (i.e. changes in slope).
summary(Diversity_mod1)

# do pairwise comparisons of intercepts
emm_intercept <- emmeans(Diversity_mod1, "Stress", at = list(log10_Diversity_centered= 0))
pairs(emm_intercept, adjust = "none" ) 

# compare slopes between treatments
slope_comp <- emtrends(Diversity_mod1, pairwise ~ Stress, var = "log10_Diversity_centered",
                       adjust = NULL
)

slope_comp$contrasts

### Plot the model
newX <- expand.grid(
  log10_Diversity_centered= seq(-0.451545 ,0.451545, length = 100),
  Stress= unique(All_data_noEvo$Stress))


# new Y's - one for average (fixed) and one for each OTU (re.form = ~....)
fixed_pred <- predict(Diversity_mod1, newdata = newX, re.form = NA)

# housekeeping - create new dataframe with predicted values
pd <- data.frame(newX, fixed_pred)

pd_plot <- pd %>%
  mutate(Stress = case_when(
    Stress == "pH_Sal" ~ "pH x Salinity",
    Stress == "pH_Temp" ~ "pH x Temperature",
    Stress == "Temp" ~ "Temperature",
    Stress =="Sal" ~ "Salinity",
    Stress =="pH" ~ "pH",
    Stress == "Control" ~ "Control",
    Stress =="Sal_Temp" ~ "Salinity x Temperature",
    Stress =="pH_Sal_Temp"~ "pH x Salinity x Temperature",
    TRUE ~ Stress 
  ))

meanBiomass_plot <- mean_Biomass %>%
  mutate(Stress = case_when(
    Stress == "pH_Sal" ~ "pH x Salinity",
    Stress == "pH_Temp" ~ "pH x Temperature",
    Stress == "Temp" ~ "Temperature",
    Stress =="Sal" ~ "Salinity",
    Stress =="pH" ~ "pH",
    Stress == "Control" ~ "Control",
    Stress =="Sal_Temp" ~ "Salinity x Temperature",
    Stress =="pH_Sal_Temp"~ "pH x Salinity x Temperature",
    TRUE ~ Stress 
  ))

All_data_noEvo_plot <- All_data_noEvo %>%
  mutate(Stress = case_when(
    Stress == "pH_Sal" ~ "pH x Salinity",
    Stress == "pH_Temp" ~ "pH x Temperature",
    Stress == "Temp" ~ "Temperature",
    Stress =="Sal" ~ "Salinity",
    Stress =="pH" ~ "pH",
    Stress == "Control" ~ "Control",
    Stress =="Sal_Temp" ~ "Salinity x Temperature",
    Stress =="pH_Sal_Temp"~ "pH x Salinity x Temperature",
    TRUE ~ Stress 
  ))

# set the theme
theme_mine <- function(base_size = 12, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16),
      axis.text.x = element_text(size=14),
      axis.text.y = element_text(size=16,hjust=1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size = 16, margin = margin(t = 20)),
      axis.title.y= element_text(size=16,angle=90),
      panel.background = element_blank(), 
      panel.border =element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}
# Plot the results
BEF_Plot_all <- ggplot(pd_plot, aes(x = log10_Diversity_centered, y = fixed_pred))+
  geom_line(linewidth = 1.5, col = "red")+
  geom_point(data = meanBiomass_plot, aes(x =log10_Diversity_centered, y = log_mean_biomass), size = 3)+
  geom_point(data = All_data_noEvo_plot, aes(x =log10_Diversity_centered, y = log10_OD), alpha = 0.2)+
  theme_mine() +
  facet_wrap(~fct_relevel(Stress, "Control", "pH", "Salinity", "Temperature", "pH x Salinity", "pH x Temperature", "Salinity x Temperature", "pH x Salinity x Temperature"), 
             nrow = 1, labeller =label_wrap_gen(width = 12)) +
  labs(x = expression("Centered Species Richness: " ~ log[10] *"(S) -" ~ log[10] * "(S"[c] *")")) +
  labs(y = expression("Community Functioing: " ~ log[10] *"(Y(OD"[600]*"))"))


BEF_Plot_all

#save as set size

png(filename="BEF_lm_plot.png", width=1800, height=500, res = 100)
BEF_Plot_all
dev.off()


#############################################
# Community composition analysis (clean)
# - Ordered pairwise heatmap axes
# - Shannon box with significance brackets
# - All plots: TIFF 600 dpi
#############################################

# ============ Setup ============
pkgs <- c(
  "magrittr","dplyr","tidyr","tibble","ggplot2","ggpubr","scales","ggrepel","patchwork",
  "ape","phyloseq","vegan","pairwiseAdonis","ggtree","stringr","picante","phangorn","microeco","reshape2"
)
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!"pairwiseAdonis" %in% installed.packages()[,1]) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}

# Libraries
library(magrittr); library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(ggpubr); library(scales); library(ggrepel); library(patchwork)
library(ape); library(phyloseq); library(vegan); library(pairwiseAdonis)
library(ggtree); library(stringr); library(picante); library(phangorn); library(microeco); library(reshape2)

set.seed(123)

# ============ Data I/O ============
# Expect these files in working directory:
# - ASV.tax.tsv
# - ASV.counts.tsv
# - Sequencing_metadata2.csv
# - ASV1-13_tree.txt

taxa_table <- read.csv("ASV.tax.tsv", sep = "\t")
taxa_table <- taxa_table %>%
  mutate(Genus = ifelse(row.names(.) == "ASV_10", "Pseudomonas2", Genus)) %>%
  microeco::tidy_taxonomy() %>%
  as.matrix()

OTU_data <- read.csv("ASV.counts.tsv", sep = "\t")
OTU_data <- as.matrix(t(OTU_data)) # samples as rows

sampleData <- read.csv("Sequencing_metadata2.csv") %>%
  mutate(
    pH       = ifelse(Stress %in% c("pH","pHSal","pHTemp","pHSalTemp"), 5.5, 7.2),
    Salinity = ifelse(Stress %in% c("Sal","pHSal","SalTemp","pHSalTemp"), 20, 0)
  )
sampleData$Temp     <- factor(sampleData$Temp)
sampleData$Salinity <- factor(sampleData$Salinity)
sampleData$pH       <- factor(sampleData$pH)
rownames(sampleData) <- sampleData$SampleID
sampleData <- dplyr::select(sampleData, -SampleID)

tree_data <- ape::read.tree("ASV1-13_tree.txt")

# ============ Phyloseq object & filtering ============
ps <- phyloseq(
  otu_table(OTU_data, taxa_are_rows = TRUE),
  sample_data(sampleData),
  tax_table(taxa_table),
  phy_tree(tree_data)
)

# Relative abundance + filter mean abundance > 0.5%
psr        <- transform_sample_counts(ps, function(x) x / sum(x))
ps_filter1 <- filter_taxa(psr, function(x) mean(x) > 0.005, prune = TRUE)

# Optional: merge two Pseudomonas if present at positions 8:9
if (ntaxa(ps_filter1) >= 9) {
  ps_filter1 <- merge_taxa(ps_filter1, taxa_names(ps_filter1)[8:9], 2)
}

# Keep only Diversity == "Eight" and non-evo
ps_filter1 <- subset_samples(ps_filter1, Diversity == "Eight" & Evo_status != "Evo")

# ============ Matrices & metadata (aligned) ============
ASV_t        <- as.data.frame(t(otu_table(ps_filter1)))                  # samples × ASVs
env_option_A <- as(phyloseq::sample_data(ps_filter1), "data.frame")      # ensure data.frame
env_option_A$Id <- rownames(env_option_A)

common_samples     <- intersect(rownames(ASV_t), rownames(env_option_A))
ASV_t_clean        <- as.matrix(ASV_t[common_samples, , drop = FALSE]); storage.mode(ASV_t_clean) <- "double"
env_option_A_clean <- env_option_A[common_samples, , drop = FALSE]

# ============ Unified Stress labels & palette ============
# Short codes kept in data; "pretty" labels used for plot axes
levels_short <- c("Control","pH","Sal","Temp","pHSal","pHTemp","SalTemp","pHSalTemp")
pretty_map   <- c(
  "Control"   = "Control",
  "pH"        = "pH",
  "Sal"       = "Salinity",
  "Temp"      = "Temperature",
  "pHSal"     = "pH x Salinity",
  "pHTemp"    = "pH x Temperature",
  "SalTemp"   = "Salinity x Temperature",
  "pHSalTemp" = "pH x Salinity x Temperature"
)
pretty_levels <- unname(pretty_map[levels_short])

env_option_A_clean$Stress <- env_option_A_clean$Stress %>%
  as.character() %>%
  {replace(., . == "Salinity", "Sal")} %>%
  {replace(., . == "Temperature", "Temp")} %>%
  factor(levels = levels_short)

pal_short <- c(
  "Control"="darkgrey","pH"="#56B4E9","Sal"="#F0E442","Temp"="#FF7043",
  "pHSal"="#009E73","pHTemp"="#CC79A7","SalTemp"="chocolate4","pHSalTemp"="#3949AB"
)
pal_pretty <- setNames(unname(pal_short), pretty_levels) # same colors for pretty labels

# ============ Distance matrix ============
Distances <- vegdist(ASV_t_clean, method = "bray")

# ============ db-RDA (gradients + interactions) ============
env_dbRDA <- env_option_A_clean %>%
  mutate(
    pH_num       = as.numeric(as.character(pH)),
    Salinity_num = as.numeric(as.character(Salinity)),
    Temp_num     = as.numeric(as.character(Temp))
  )

dbRDA_grad <- capscale(Distances ~ pH_num * Salinity_num * Temp_num, data = env_dbRDA)
site_scores_grad <- scores(dbRDA_grad, display = "sites") %>% as.data.frame()
site_scores_grad$Stress <- env_option_A_clean$Stress
var_expl_grad <- {
  eig <- eigenvals(dbRDA_grad, constrained = TRUE)
  round(100 * eig / sum(eig), 1)
}

dbRDA_plot_grad <- ggplot(site_scores_grad, aes(CAP1, CAP2, colour = Stress)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Stress, colour = Stress),
               type = "t", level = 0.95,
               linewidth = 0.8, alpha = 0.25) +
  labs(x = paste0("RDA1 (", var_expl_grad[1], "%)"),
       y = paste0("RDA2 (", var_expl_grad[2], "%)")) +
  scale_colour_manual(values = pal_short,
                      limits = levels_short,
                      breaks = levels_short) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.title   = element_text(size = 16, family = "Helvetica"),
    axis.text    = element_text(size = 14, family = "Helvetica"),
    legend.title = element_text(size = 16, family = "Helvetica", face = "bold"),
    legend.text  = element_text(size = 14, family = "Helvetica"),
    legend.position = "right"
  )

ggsave("dbRDA_gradients.tiff", dbRDA_plot_grad,
       width = 18, height = 15, units = "cm",
       dpi = 600, device = ragg::agg_tiff, compression = "lzw")


# ============ PERMANOVA (global + pairwise heatmap, ordered axes) ============
global_perm <- adonis2(Distances ~ Stress, data = env_option_A_clean, permutations = 999)
cat("Global PERMANOVA: R² =", round(global_perm$R2[1], 2),
    ", p =", signif(global_perm$`Pr(>F)`[1], 3), "\n")

pw <- pairwise.adonis(x = Distances, factors = env_option_A_clean$Stress, perm = 999) %>%
  as.data.frame()

# Parse groups (short codes)
pw$Group1 <- sub(" vs .*", "", pw$pairs)
pw$Group2 <- sub(".* vs ", "", pw$pairs)

# Determine present levels (keep specified order)
present_short <- levels_short[levels_short %in% unique(c(as.character(env_option_A_clean$Stress)))]
# Initialize full R2 matrix in that order
mat <- matrix(NA_real_, nrow = length(present_short), ncol = length(present_short),
              dimnames = list(present_short, present_short))
# Fill upper triangle from results
for (i in seq_len(nrow(pw))) {
  g1 <- pw$Group1[i]; g2 <- pw$Group2[i]
  if (g1 %in% present_short && g2 %in% present_short) {
    mat[g1, g2] <- pw$R2[i]
  }
}
# Symmetrize
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
# Long format
mat_long <- reshape2::melt(mat, varnames = c("Group1","Group2"), value.name = "R2", na.rm = TRUE)

# add adjusted p and stars (robust to column name variants)
pcol <- intersect(c("p.adjusted","p.adjustment","p.adj"), names(pw))
if (length(pcol) == 0) stop("Could not find adjusted p-value column in pairwise.adonis output.")
sig_map <- pw[, c("Group1","Group2", pcol[1])]
names(sig_map)[3] <- "p.adj"
sig_map$stars <- cut(sig_map$p.adj,
                     breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                     labels = c("***","**","*",""))

mat_long <- mat_long %>% left_join(sig_map, by = c("Group1","Group2"))

# Recode axes to pretty labels in the same order
pretty_present <- unname(pretty_map[present_short])
mat_long <- mat_long %>%
  mutate(
    Group1 = recode(Group1, !!!as.list(pretty_map)),
    Group2 = recode(Group2, !!!as.list(pretty_map)),
    Group1 = factor(Group1, levels = pretty_present),
    Group2 = factor(Group2, levels = pretty_present)
  )

heatmap_plot <- ggplot(mat_long, aes(x = Group1, y = Group2, fill = R2)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkred", na.value = "grey90") +
  geom_text(aes(label = stars), vjust = 0.5, size = 5, color = "black", na.rm = TRUE) +
  labs(x = NULL, y = NULL, fill = expression(R^2)) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x   = element_text(size = 14, angle = 45, hjust = 1, vjust = 1, family = "Helvetica"),
    axis.text.y   = element_text(size = 14, family = "Helvetica"),  # upright
    axis.title    = element_text(size = 16, family = "Helvetica"),
    legend.title  = element_text(size = 16, face = "bold", family = "Helvetica"),
    legend.text   = element_text(size = 14, family = "Helvetica")
  )

ggsave("PERMANOVA_pairwise_heatmap.tiff", heatmap_plot,
       width = 18, height = 15, units = "cm",
       dpi = 600, device = ragg::agg_tiff, compression = "lzw")

# Show all pairwise PERMANOVA results as a data.frame
print(pw)

readr::write_csv(tibble::rownames_to_column(as.data.frame(pw), "contrast"),
                 "Supp_Table_PERMANOVA_pairwise.csv")


# ============================
# Shannon diversity (Nature-ready, robust Id handling)
# ============================

# Core tidy + plotting
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(scales)

# Stats + tables
library(vegan)
library(broom)
library(car)          # Anova(), leveneTest()  (NOTE: masks dplyr::recode if unqualified)
library(effectsize)   # eta_squared(), omega_squared()
library(emmeans)      # emmeans(), pairs(), eff_size()
library(readr)

# ---------- Preconditions / sanity checks ----------
# env_option_A_clean: metadata with Id + Stress
# ASV_t_clean:       OTU/ASV matrix aligned to samples (columns)
# levels_short:      factor levels in canonical short code order for Stress
# pretty_map:        named character vector mapping short -> pretty labels
# pretty_levels:     ordered vector of pretty labels for plotting
# pal_pretty:        named palette for pretty labels

# Ensure 'Id' exists exactly once in env_option_A_clean
if (!"Id" %in% names(env_option_A_clean)) {
  env_option_A_clean <- env_option_A_clean %>% tibble::rownames_to_column("Id")
}

# Basic validation of pretty_map / pretty_levels
if (!exists("pretty_map") || !is.character(pretty_map) || is.null(names(pretty_map))) {
  stop("`pretty_map` must be a *named* character vector mapping short levels to pretty labels.")
}
if (!exists("pretty_levels") || !is.character(pretty_levels)) {
  stop("`pretty_levels` must be a character vector giving the desired order of pretty labels.")
}

# Helper: safely pull first existing column from a data.frame
grab_col <- function(df, candidates, default = NA_real_) {
  for (nm in candidates) {
    if (nm %in% names(df)) return(df[[nm]])
  }
  rep(default, nrow(df))
}

# ---------- 1) Shannon on aligned matrix ----------
Data_richness <- vegan::diversity(ASV_t_clean, index = "shannon")

# Tidy richness with sample IDs
Data_richness_df <- tibble(
  Id = names(Data_richness),
  shan_div = as.numeric(Data_richness)
)

# ---------- 2) Join to metadata and add pretty labels ----------
div_df <- env_option_A_clean %>%
  select(Id, Stress) %>%
  left_join(Data_richness_df, by = "Id") %>%
  mutate(
    Stress        = factor(as.character(Stress), levels = levels_short),     # ensure canonical order
    # Qualify dplyr::recode to avoid clash with car::recode; keep unmapped as default
    Stress_pretty = dplyr::recode(as.character(Stress), !!!pretty_map, .default = as.character(Stress)),
    Stress_pretty = factor(Stress_pretty, levels = pretty_levels)
  )

# Optional sanity check: warn if any levels weren’t mapped
unmapped <- setdiff(levels(div_df$Stress), names(pretty_map))
if (length(unmapped) > 0) {
  warning("These Stress levels were not found in pretty_map and will be left as-is: ",
          paste(unmapped, collapse = ", "))
}

# ---------- 3) ANOVA + Tukey (original objects for plotting) ----------
Shannon_anova <- aov(shan_div ~ Stress, data = div_df)

tukey_df <- TukeyHSD(Shannon_anova, conf.level = 0.95)[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("pairwise")
names(tukey_df) <- make.names(names(tukey_df))  # "p adj" -> "p.adj"

# ---------- Plot (boxplot) ----------
plot_shandiv_box <- ggplot(div_df, aes(x = Stress_pretty, y = shan_div, fill = Stress_pretty)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.6) +
  scale_fill_manual(values = pal_pretty, limits = pretty_levels) +
  scale_x_discrete(labels = scales::label_wrap(12)) +
  labs(x = NULL, y = "Shannon diversity (H′)") +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x  = element_text(angle = 50, hjust = 1, vjust = 1, size = 14, family = "Helvetica"),
    axis.text.y  = element_text(size = 14, family = "Helvetica"),
    axis.title.y = element_text(size = 16, family = "Helvetica"),
    legend.position = "none",
    plot.margin = margin(8, 12, 8, 8)
  )

# ---------- 4) Nature-ready tables & diagnostics (keeps your plot code unchanged) ----------
# Fit as lm for Type II SS (robust to imbalance when no interactions are modeled)
fit_shan <- lm(shan_div ~ Stress, data = div_df)

# Assumptions (report in Supplementary Methods)
assump_shapiro <- shapiro.test(residuals(fit_shan))                   # residual normality
assump_levene  <- car::leveneTest(shan_div ~ Stress, div_df)          # homogeneity of var
welch_alt      <- oneway.test(shan_div ~ Stress, div_df, var.equal = FALSE)  # sensitivity

# Type II ANOVA table
anova_type2 <- car::Anova(fit_shan, type = 2) %>% broom::tidy()

# Effect sizes with CIs (function returns different column names depending on context; handle both)
eta2_raw <- effectsize::eta_squared(fit_shan, partial = FALSE, ci = 0.95) %>% as.data.frame()
omega_raw <- effectsize::omega_squared(fit_shan, partial = FALSE) %>% as.data.frame()

eta2_tbl <- tibble(
  term   = grab_col(eta2_raw, c("Effect", "Term", "Parameter"), default = NA_character_),
  eta2   = grab_col(eta2_raw, c("Eta2_partial", "Eta2")),
  ci_low = grab_col(eta2_raw, c("CI_low", "CI_low_partial", "CI_low_Eta2")),
  ci_high= grab_col(eta2_raw, c("CI_high", "CI_high_partial", "CI_high_Eta2"))
) %>% distinct()

omega2_tbl <- tibble(
  term   = grab_col(omega_raw, c("Effect", "Term", "Parameter"), default = NA_character_),
  omega2 = grab_col(omega_raw, c("Omega2_partial", "Omega2"))
) %>% distinct()

S1a <- anova_type2 %>%
  filter(term == "Stress") %>%
  transmute(
    term,
    df     = as.integer(df),
    ss     = sumsq,
    ms     = sumsq/df,
    F      = statistic,
    p      = p.value,
    N      = nrow(div_df),
    k      = dplyr::n_distinct(div_df$Stress)
  ) %>%
  left_join(eta2_tbl,  by = "term") %>%
  left_join(omega2_tbl, by = "term")

S1a_print <- S1a %>%
  mutate(across(c(ss, ms, F, eta2, omega2, ci_low, ci_high), ~round(., 3)),
         p = signif(p, 3))

# --- Post-hoc contrasts (Tukey) with Hedges' g and 95% CIs ---
emm <- emmeans::emmeans(fit_shan, ~ Stress)

# Pairwise differences
pw_contr <- contrast(emm, method = "pairwise")

# Tukey-adjusted mean-difference table (ensure p-values are returned)
tuk_raw <- summary(pw_contr, infer = c(FALSE, TRUE), adjust = "tukey") %>% as.data.frame()

# Pick the stat column (t.ratio for LM; z.ratio for some GLMs)
stat_col <- if ("t.ratio" %in% names(tuk_raw)) "t.ratio" else
  if ("z.ratio" %in% names(tuk_raw)) "z.ratio" else NA_character_

tuk <- tuk_raw %>%
  dplyr::transmute(
    contrast = .data$contrast,
    estimate = .data$estimate,   # mean difference
    SE       = .data$SE,
    df       = if ("df" %in% names(tuk_raw)) .data$df else NA_real_,
    stat     = if (!is.na(stat_col)) .data[[stat_col]] else NA_real_,
    p.adj    = if ("p.value" %in% names(tuk_raw)) .data$p.value else NA_real_
  )

# --- Hedges' g computed manually from mean diffs ---
sigma_hat <- sigma(fit_shan)
df_resid  <- df.residual(fit_shan)
J_corr    <- 1 - 3/(4*df_resid - 1)              # small-sample correction

S1b <- tuk %>%
  dplyr::mutate(
    hedges_g = J_corr * (estimate / sigma_hat),
    SE_g     = J_corr * (SE / sigma_hat),
    t_crit   = ifelse(is.finite(df), qt(0.975, df), NA_real_),
    g_low    = hedges_g - t_crit * SE_g,
    g_high   = hedges_g + t_crit * SE_g
  ) %>%
  dplyr::select(contrast, estimate, SE, df, stat, p.adj, hedges_g, g_low, g_high) %>%
  dplyr::arrange(p.adj)

# Write supplementary tables
readr::write_csv(S1a_print, "Shannon_ANOVA_effectsizes.csv")
readr::write_csv(S1b,       "Tukey_posthoc_effectsizes.csv")


# ---------- 5) Significance brackets: Control vs each; only significant (your original logic) ----------
sig_df <- tukey_df %>%
  tidyr::separate(pairwise, into = c("g1","g2"), sep = "-") %>%
  dplyr::mutate(
    g1 = dplyr::recode(g1, !!!pretty_map, .default = g1),
    g2 = dplyr::recode(g2, !!!pretty_map, .default = g2)
  ) %>%
  dplyr::filter(g1 == "Control" | g2 == "Control") %>%
  dplyr::arrange(p.adj) %>%
  dplyr::filter(p.adj < 0.05)

plot_shandiv_box_sig <- plot_shandiv_box
if (nrow(sig_df) > 0) {
  y_range <- range(div_df$shan_div, na.rm = TRUE)
  step <- 0.06 * diff(y_range)
  base <- max(div_df$shan_div, na.rm = TRUE) + step
  
  sig_df <- sig_df %>%
    dplyr::transmute(
      group1 = ifelse(g1 == "Control", g1, g2),
      group2 = ifelse(g1 == "Control", g2, g1),
      p.adj,
      y.position = base + dplyr::row_number() * step,
      p.signif = dplyr::case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        TRUE          ~ "*"
      )
    )
  
  plot_shandiv_box_sig <- plot_shandiv_box +
    ggpubr::stat_pvalue_manual(
      sig_df,
      label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.3,
      size = 3
    ) +
    expand_limits(y = max(sig_df$y.position) + step)
}

# ---------- 6) Save (TIFF, 600 dpi) ----------
ggsave("Shannon_diversity_boxplot_with_signif.tiff", plot_shandiv_box_sig,
       width = 18, height = 15, units = "cm", dpi = 600, compression = "lzw")

# ---------- Console summary ----------
cat("Saved:\n",
    "- dbRDA_gradients.tiff\n",
    "- PERMANOVA_pairwise_heatmap.tiff\n",
    "- Shannon_diversity_boxplot_with_signif.tiff\n",
    "- Combined_Figure_3panels.tiff (optional)\n")

cat("\nANOVA (Type II) summary for manuscript:\n")
cat(sprintf("F(%d, %d) = %.2f, p = %s; partial eta^2 = %.3f [%0.3f, %0.3f]; partial omega^2 = %.3f.\n",
            S1a$df, df.residual(fit_shan), S1a$F, signif(S1a$p, 3),
            S1a$eta2, S1a$ci_low, S1a$ci_high, S1a$omega2))

cat("\nAssumptions (report in Supplementary Methods):\n")
cat(sprintf("- Shapiro-Wilk on residuals: W=%.3f, p=%s\n", assump_shapiro$statistic, signif(assump_shapiro$p.value,3)))
cat(sprintf("- Levene's test: F=%.2f, p=%s\n", assump_levene$`F value`[1], signif(assump_levene$`Pr(>F)`[1],3)))
cat(sprintf("- Welch ANOVA (sensitivity): F=%.2f, p=%s\n", welch_alt$statistic, signif(welch_alt$p.value,3)))



# ============================
# Mantel test (taxa-level)
# ============================

# 1) Root the tree at midpoint (class: "phylo")
tree_rooted <- phangorn::midpoint(phy_tree(ps_filter1))

# 2) Ecological distances between taxa (ASVs): Bray on transposed ASV table
dist_ecol_m <- as.matrix(vegan::vegdist(t(ASV_t_clean), method = "bray"))

# 3) Phylogenetic distances between taxa (patristic)
#    Use cophenetic.phylo() explicitly (some ape builds hide the generic)
if ("cophenetic.phylo" %in% getNamespaceExports("ape")) {
  dist_phylo_m <- ape::cophenetic.phylo(tree_rooted)
} else if ("cophenetic" %in% getNamespaceExports("ape")) {
  dist_phylo_m <- ape::cophenetic(tree_rooted)
} else {
  stop("Your 'ape' installation does not expose cophenetic methods. Try reinstalling/updating 'ape'.")
}

# 4) Align the two distance matrices to common ASV IDs
common_taxa <- intersect(rownames(dist_ecol_m), rownames(dist_phylo_m))
dist_ecol_taxa <- as.dist(dist_ecol_m[common_taxa, common_taxa, drop = FALSE])
dist_phylo     <- as.dist(dist_phylo_m[common_taxa, common_taxa, drop = FALSE])

# 5) Mantel test
set.seed(123)
mantel_res <- vegan::mantel(dist_ecol_taxa, dist_phylo,
                            method = "spearman", permutations = 9999)
print(mantel_res)

# 6) Mantel scatter (Nature-style)
mantel_df <- data.frame(
  Phylogenetic = as.vector(dist_phylo),
  Ecological   = as.vector(dist_ecol_taxa)
)

# Mantel stats (numeric)
mantel_r_num <- as.numeric(mantel_res$statistic)
mantel_p_num <- as.numeric(mantel_res$signif)

# Safe strings for display (always non-empty)
r_text <- formatC(mantel_r_num, format = "f", digits = 3)
p_text <- if (is.finite(mantel_p_num)) {
  formatC(mantel_p_num, format = "e", digits = 2)   # scientific form
} else {
  "NA"
}

# Plain label string (no plotmath)
label_str <- paste0("r = ", r_text, ", p = ", p_text)

mantel_plot <- ggplot(mantel_df, aes(x = Phylogenetic, y = Ecological)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 0.8) +
  theme_classic(base_size = 14, base_family = "Helvetica") +   # axis text = 14
  labs(
    x = "Phylogenetic distance",
    y = "Ecological distance (Bray–Curtis)"
  ) +
  theme(
    axis.title = element_text(size = 16, family = "Helvetica") # axis titles = 16
  ) +
  annotate(
    "text",
    x = quantile(mantel_df$Phylogenetic, 0.05, na.rm = TRUE),
    y = quantile(mantel_df$Ecological,   0.95, na.rm = TRUE),
    label = label_str,
    hjust = 0, vjust = 1,
    size = 5,                  # ≈14 pt in ggplot units
    family = "Helvetica"
  )



ggsave("Mantel_taxa_plot.tiff", mantel_plot,
       width = 18, height = 15, units = "cm", dpi = 600, compression = "lzw")


# =========================================================
# ISA — Top-1 indicator per stress (one-vs-rest) + Nature-style plot
# Single colourbar (Indicator value), constant point size
# =========================================================
suppressPackageStartupMessages({
  library(indicspecies)
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(ggplot2); library(scales)
  library(permute)
})

set.seed(123)

# ---- Required inputs in memory ----
# ASV_t_clean        # samples x ASVs (rows = samples, cols = ASVs)
# env_option_A_clean # data.frame with column 'Stress' (one per sample, same row order as ASV_t_clean)
# taxa_table         # taxonomy; rownames = ASV, has a 'Genus' column (fallback to ASV if missing)

grp <- env_option_A_clean$Stress
stopifnot(nrow(ASV_t_clean) == length(grp))

# ---- Pretty labels / order (fallback if you don’t already have these) ----
if (!exists("pretty_map")) {
  levels_short <- c("Control","pH","Sal","Temp","pHSal","pHTemp","SalTemp","pHSalTemp")
  pretty_map <- c(
    "Control"="Control",
    "pH"="pH",
    "Sal"="Salinity",
    "Temp"="Temperature",
    "pHSal"="pH × Salinity",
    "pHTemp"="pH × Temperature",
    "SalTemp"="Salinity × Temperature",
    "pHSalTemp"="pH × Salinity × Temperature"
  )
}
short_levels   <- names(pretty_map)
present_short  <- short_levels[short_levels %in% unique(as.character(grp))]
present_pretty <- unname(pretty_map[present_short])

# ---- Taxonomy table ----
tax_df <- as.data.frame(taxa_table) %>% tibble::rownames_to_column("ASV")

# ---- One-vs-rest ISA (no permutations) ----
TopN       <- 1
assoc_func <- "IndVal.g"

panel_list <- lapply(present_short, function(g) {
  grp_bin <- factor(ifelse(as.character(grp) == g, g, "Other"),
                    levels = c("Other", g))
  
  fit <- indicspecies::multipatt(
    ASV_t_clean, grp_bin,
    func = assoc_func, max.order = 1, duleg = FALSE,
    control = permute::how(nperm = 0)
  )
  
  sig <- as.data.frame(fit$sign); sig$ASV <- rownames(sig)
  mem_col <- paste0("s.", levels(grp_bin)[2])  # "s.<g>"
  
  sig %>%
    transmute(ASV, stat = as.numeric(stat), mem = .data[[mem_col]]) %>%
    filter(is.finite(stat), stat > 0, mem == 1) %>%
    select(-mem) %>%
    left_join(tax_df, by = "ASV") %>%
    mutate(
      Genus = stringr::str_remove(coalesce(Genus, ""), "^g__"),
      Genus = ifelse(Genus == "" | is.na(Genus), ASV, Genus),
      Genus = dplyr::case_when(
        Genus == "Pseudomonas"  ~ "Pseudomonas.A",
        Genus == "Pseudomonas2" ~ "Pseudomonas.F",
        TRUE ~ Genus
      ),
      Group_short  = g,
      Group_pretty = pretty_map[g]
    ) %>%
    arrange(desc(stat))
})

# ---- Keep both full stats and TopN ----
all_isa_stats <- bind_rows(panel_list) %>%
  mutate(Group_pretty = factor(Group_pretty, levels = present_pretty))

topN_by_group <- all_isa_stats %>%
  group_by(Group_short) %>%
  slice_head(n = TopN) %>%
  ungroup()

# (optional) write full stats to file for inspection
write.csv(all_isa_stats, "ISA_stats_full.csv", row.names = FALSE)

# =========================================================
# Nature-style plot: constant point size, single colourbar
# =========================================================

# Colourbar breaks & labels (nice 3 ticks)
rng   <- range(topN_by_group$stat, na.rm = TRUE)
brks0 <- pretty(rng, n = 3)
brks  <- if (length(brks0) >= 3) brks0[c(1, floor(length(brks0)/2), length(brks0))] else rng
labf  <- scales::number_format(accuracy = 0.01)

# Blue gradient (swap to viridis if you prefer)
col_vec <- c("#e8f1fb", "#6aa0d8", "#08306b")

POINT_SIZE <- 5.5  # constant, readable size

p_isa <- ggplot(topN_by_group, aes(x = Group_pretty, y = Genus)) +
  geom_point(aes(colour = stat), size = POINT_SIZE, shape = 16, alpha = 0.9) +
  scale_colour_gradientn(
    colours = col_vec,
    limits  = rng,
    breaks  = brks,
    labels  = labf,
    name    = "Indicator value"
  ) +
  guides(
    colour = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(80, "pt"), barheight = unit(6, "pt"),
      ticks.colour = "grey20"
    )
  ) +
  scale_x_discrete(labels = scales::label_wrap(14)) +
  labs(x = "Stress condition", y = "Indicator taxon (Genus)") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position   = "top",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key        = element_rect(fill = "white", colour = NA),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 12),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y       = element_text(size = 14, face = "italic"),
    axis.title.x      = element_text(size = 16),
    axis.title.y      = element_text(size = 16),
    axis.line         = element_line(linewidth = 0.4),
    axis.ticks.length = unit(3, "pt"),
    plot.margin       = margin(8, 12, 8, 8)
  )

ggsave("ISA_Top1_1vsRest_singlepanel.tiff", p_isa,
       width = 18, height = 15, units = "cm",
       dpi = 600, compression = "lzw")


# ============================================================
# log10(Community yield) ~ variance in response traits (growth)
# Nature/Science-style plots (Helvetica, no title), S = 2 removed
# Needs in working dir:
# - CommunityOD_All_Jan23_edit_zeros.csv
# - CommunityAbundance_Growth_Data.csv
# - LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv
# Optional:
# - LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv
# - LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv
# Outputs: PNG & TIFF (unweighted + weighted) + side-by-side panel
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c("tidyverse","broom","ggpubr","rlang","readr","scales")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

# ---------- Fonts (Windows fallback so Helvetica -> Arial) ----------
if (.Platform$OS.type == "windows") {
  if (!"Helvetica" %in% names(grDevices::windowsFonts())) {
    grDevices::windowsFonts(Helvetica = grDevices::windowsFont("Arial"))
  }
}

# ---------- Helpers ----------
to_short_stress <- function(x) {
  y <- gsub("_", "", as.character(x), fixed = TRUE)
  y[y == "Salinity"]    <- "Sal"
  y[y == "Temperature"] <- "Temp"
  y
}
levels_short <- c("Control","pH","Sal","Temp","pHSal","pHTemp","SalTemp","pHSalTemp")

wvar_unbiased <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  sw <- sum(w); if (sw <= 0) return(NA_real_)
  mu  <- sum(w * x) / sw
  num <- sum(w * (x - mu)^2)
  den <- sw - (sum(w^2) / sw)
  if (den <= 0) return(NA_real_)
  num / den
}
read_repair <- function(path) read.csv(path, check.names = FALSE) |> tibble::as_tibble(.name_repair = "universal")
derive_CommunityNum <- function(df, cols = c("Community","CommunityNum","Community_Number","Community.Num",
                                             "Com_Id","ComID","CommunityID")) {
  num_from_col <- function(nm) {
    if (!nm %in% names(df)) return(rep(NA_integer_, nrow(df)))
    v_chr <- as.character(df[[nm]])
    v_chr[trimws(v_chr) %in% c("", "Blank", "blank", "BLANK")] <- NA_character_
    if (is.numeric(df[[nm]])) as.integer(df[[nm]]) else suppressWarnings(as.integer(readr::parse_number(v_chr)))
  }
  cands <- lapply(cols, num_from_col)
  Reduce(dplyr::coalesce, cands) |> as.integer()
}

# ---------- 1) Growth (Id × Stress) ----------
growth_list <- list()

# pH growth (present)
d_pH <- read_repair("LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv") |>
  rename(
    Id            = any_of(c("Id","ID","id")),
    estimate      = any_of(c("estimate","Estimate","est","r","growth_rate")),
    TempLevel     = any_of(c("TempLevel","Temp","Temperature","TemperatureLevel","Temp_Level")),
    SalinityLevel = any_of(c("SalinityLevel","Salinity","Salinity_Level","Salinity.Level")),
    pH_Level      = any_of(c("pH_Level","pH","pHLevel","pH.Level"))
  ) |>
  select(Id, estimate, TempLevel, SalinityLevel, pH_Level)

d_pH2 <- d_pH |>
  mutate(Stress = case_when(
    SalinityLevel == 0  & TempLevel == 20 ~ "Control",
    SalinityLevel == 0  & TempLevel == 38 ~ "Temp",
    SalinityLevel == 20 & TempLevel == 20 ~ "Sal",
    SalinityLevel == 20 & TempLevel == 38 ~ "SalTemp"
  )) |>
  filter(pH_Level %in% c(5.5, 7.2)) |>
  mutate(
    Stress = ifelse(pH_Level == 5.5, paste0("pH", Stress), Stress),
    Stress = ifelse(Stress == "pHControl", "pH", Stress)
  ) |>
  select(Id, Stress, estimate) |>
  mutate(Stress = as.character(Stress))
growth_list[["pH"]] <- d_pH2

# Temp growth (optional)
if (file.exists("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv")) {
  d_temp <- read_repair("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv") |>
    rename(Id=any_of(c("Id","ID","id")),
           Stress=any_of(c("Stress","Treatment","Env","Condition")),
           TempLevel=any_of(c("TempLevel","Temp","Temperature","TemperatureLevel","Temp_Level")),
           estimate=any_of(c("estimate","Estimate","est","r","growth_rate"))) |>
    select(Id, Stress, TempLevel, estimate) |>
    filter(TempLevel %in% c(20, 38)) |>
    mutate(
      Stress = ifelse(TempLevel == 38, paste0(Stress,"Temp"), Stress),
      Stress = dplyr::recode(Stress, "ControlTemp"="Temp", "pH_SalTemp"="pHSalTemp", "pH_Sal"="pHSal")
    ) |>
    select(Id, Stress, estimate) |>
    mutate(Stress = as.character(Stress))
  growth_list[["Temp"]] <- d_temp
}

# Sal growth (optional)
if (file.exists("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv")) {
  d_sal <- read_repair("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv") |>
    rename(Id=any_of(c("Id","ID","id")),
           pH_Level=any_of(c("pH_Level","pH","pHLevel","pH.Level")),
           TempLevel=any_of(c("TempLevel","Temp","Temperature","TemperatureLevel","Temp_Level")),
           SalinityLevel=any_of(c("SalinityLevel","Salinity","Salinity_Level","Salinity.Level")),
           estimate=any_of(c("estimate","Estimate","est","r","growth_rate"))) |>
    select(Id, pH_Level, TempLevel, SalinityLevel, estimate) |>
    mutate(Stress = case_when(
      pH_Level == 7.2 & TempLevel == 20 ~ "Control",
      pH_Level == 7.2 & TempLevel == 38 ~ "Temp",
      pH_Level == 5.5 & TempLevel == 20 ~ "pH",
      pH_Level == 5.5 & TempLevel == 38 ~ "pHTemp"
    )) |>
    filter(SalinityLevel %in% c(0, 20)) |>
    mutate(
      Stress = ifelse(SalinityLevel == 20, paste0(Stress,"Sal"), Stress),
      Stress = dplyr::recode(Stress, "ControlSal"="Sal", "pHTempSal"="pHSalTemp", "TempSal"="SalTemp")
    ) |>
    select(Id, Stress, estimate) |>
    mutate(Stress = as.character(Stress))
  growth_list[["Sal"]] <- d_sal
}

All_growth <- bind_rows(growth_list) |> mutate(Stress = as.character(Stress))
Growth_mean <- All_growth |>
  group_by(Id, Stress) |>
  summarise(meanGrowth = mean(estimate, na.rm = TRUE), .groups = "drop")

# ---------- 2) Abundance per sample ----------
abund_raw <- read_repair("CommunityAbundance_Growth_Data.csv")
abund <- abund_raw |>
  rename(
    SampleID  = any_of(c("SampleID","sampleID","sample_id","Sample.Id")),
    Com_Id    = any_of(c("Com_Id","ComID","CommunityID","Community_Id","ComId")),
    Community = any_of(c("Community","CommunityNum","Community_Number","Community.Num")),
    Stress    = any_of(c("Stress","Treatment","Env","Condition")),
    Diversity = any_of(c("Diversity","Diversity_Level","DiversityLevel")),
    Id        = any_of(c("Id","ID","id")),
    Abundance = any_of(c("Abundance","abundance","RelAbundance","RelativeAbundance","rel_abund"))
  ) |>
  mutate(
    CommunityNum = derive_CommunityNum(cur_data()),
    Stress = to_short_stress(Stress)
  ) |>
  select(SampleID, Com_Id, CommunityNum, Stress, Diversity, Id, Abundance) |>
  filter(!is.na(CommunityNum), !is.na(Stress), !is.na(Id))

abund_g <- abund |> left_join(Growth_mean, by = c("Id","Stress"))

# ---------- 3) Trait variance at Community × Stress ----------
# Unweighted
trait_unw <- abund_g |>
  group_by(CommunityNum, Stress) |>
  summarise(
    n_taxa    = n_distinct(Id[is.finite(meanGrowth)]),
    Trait_var = var(meanGrowth, na.rm = TRUE),
    .groups = "drop"
  )

# Weighted per sample -> averaged to Community×Stress
trait_w_sample <- abund_g |>
  group_by(SampleID, CommunityNum, Stress) |>
  summarise(
    n_taxa_w = n_distinct(Id[is.finite(meanGrowth) & is.finite(Abundance) & Abundance > 0]),
    wvar     = wvar_unbiased(meanGrowth, Abundance),
    .groups = "drop"
  )
trait_w <- trait_w_sample |>
  group_by(CommunityNum, Stress) |>
  summarise(Trait_var_w = mean(wvar, na.rm = TRUE), .groups = "drop")

suppressPackageStartupMessages({
  library(dplyr); library(broom); library(emmeans); library(readr)
})

# ---------- Build analysis table (weighted) ----------
dw <- yield_cs %>%
  left_join(trait_w, by = c("CommunityNum","Stress")) %>%
  mutate(
    logVar   = log10(pmax(Trait_var_w, .Machine$double.xmin)),
    logYield = log10_Yield,
    Richness_f = factor(Richness, levels = sort(unique(Richness)))
  ) %>%
  filter(is.finite(logVar), is.finite(logYield))

# ---------- ANCOVA ----------
fit_ancova <- lm(logYield ~ logVar * Stress + Richness_f, data = dw)

ancova_tab <- anova(fit_ancova)
int_rowname <- rownames(ancova_tab)[grepl("(logVar:Stress|Stress:logVar)", rownames(ancova_tab))][1]
ancova_stats <- tibble(
  term   = int_rowname,
  df1    = ancova_tab[[ "Df"       ]][match(int_rowname, rownames(ancova_tab))],
  df2    = ancova_tab[[ "Df"       ]][nrow(ancova_tab)],   # residual df
  F      = ancova_tab[[ "F value"  ]][match(int_rowname, rownames(ancova_tab))],
  P      = ancova_tab[[ "Pr(>F)"   ]][match(int_rowname, rownames(ancova_tab))],
  R2_adj = summary(fit_ancova)$adj.r.squared
)

# ---------- Per-stress slopes (±95% CI, P) ----------
varname <- "logVar"  # the predictor you trended over
tr <- emtrends(fit_ancova, ~ Stress, var = varname)
tr_sum <- summary(tr, infer = c(TRUE, TRUE))    # adds 95% CI and P
tr_df  <- as.data.frame(tr_sum)

# emtrends names the slope column "<var>.trend" (e.g., "logVar.trend")
trend_col <- grep("\\.trend$", names(tr_df), value = TRUE)
stopifnot(length(trend_col) == 1)

slopes_emm <- tr_df %>%
  transmute(
    Stress,
    slope   = .data[[trend_col]],
    SE      = SE,
    df      = df,
    t       = t.ratio,
    P       = p.value,
    CI_low  = lower.CL,
    CI_high = upper.CL
  )

# ---------- Add n per stress and per-stress R² (from separate fits) ----------
n_by_stress <- dw %>% count(Stress, name = "n")

r2_by_stress <- dw %>%
  group_by(Stress) %>%
  do({
    m  <- lm(logYield ~ logVar, data = .)
    gl <- broom::glance(m)
    tibble(R2 = gl$r.squared, R2_adj = gl$adj.r.squared)
  }) %>% ungroup()

slopes_final <- slopes_emm %>%
  left_join(n_by_stress,  by = "Stress") %>%
  left_join(r2_by_stress, by = "Stress") %>%
  arrange(Stress)

# ---------- Pairwise differences between slopes (optional) ----------
ctr <- summary(pairs(tr), infer = c(TRUE, TRUE))
slope_contrasts <- as_tibble(ctr) %>%
  rename(contrast = contrast,
         estimate = estimate,
         SE = SE, df = df,
         t = t.ratio, P = p.value,
         CI_low = lower.CL, CI_high = upper.CL)

# ---------- Save outputs ----------
write_csv(slopes_final,    "TraitVariance_WEIGHTED_slopes_by_stress.csv")
write_csv(ancova_stats,    "TraitVariance_WEIGHTED_ANCOVA_interaction.csv")
write_csv(slope_contrasts, "TraitVariance_WEIGHTED_slope_contrasts.csv")

# Inspect in console
print(slopes_final)
print(ancova_stats)


# ---------- 4) Community yields ----------
com_raw <- read_repair("CommunityOD_All_Jan23_edit_zeros.csv")
yield_cs <- com_raw |>
  rename(
    t               = any_of(c("t","time","Timepoint")),
    Evo_Treatment   = any_of(c("Evo_Treatment","Evo_status","Evo","Evolution")),
    Stress          = any_of(c("Stress","Treatment","Env","Condition")),
    Community       = any_of(c("Community","CommunityNum","Community_Number","Community.Num")),
    Com_Id          = any_of(c("Com_Id","ComID","CommunityID","Community_Id","ComId")),
    Diversity_Level = any_of(c("Diversity_Level","Diversity","Richness")),
    OD              = any_of(c("OD","OD600","OD_600"))
  ) |>
  mutate(
    Stress = to_short_stress(Stress),
    CommunityNum = derive_CommunityNum(cur_data()),
    Richness = as.integer(Diversity_Level)
  ) |>
  filter(t == 2, Evo_Treatment != "Evo", !is.na(CommunityNum)) |>
  group_by(CommunityNum, Stress) |>
  summarise(
    Yield    = mean(OD, na.rm = TRUE),
    Richness = dplyr::first(na.omit(Richness)),
    .groups = "drop"
  ) |>
  mutate(log10_Yield = log10(pmax(Yield, .Machine$double.xmin)))

# ---------- 5) Combine (now including 2-taxa) ----------
df <- yield_cs |>
  left_join(trait_unw, by = c("CommunityNum","Stress")) |>
  left_join(trait_w,   by = c("CommunityNum","Stress")) |>
  mutate(
    Stress     = factor(Stress, levels = levels_short),
    Richness   = as.integer(Richness),
    Richness_f = factor(Richness,
                        levels = c(2,4,8),
                        labels = c("Two","Four","Eight"))
  )

# ---------- 6) Nature/Science-style plots ----------
dir.create("plots", showWarnings = FALSE)

pal <- c("Control"="grey30","pH"="#56B4E9","Sal"="#F0E442","Temp"="#FF7043",
         "pHSal"="#009E73","pHTemp"="#CC79A7","SalTemp"="chocolate4","pHSalTemp"="#3949AB")

# Updated shapes to include Two-taxa communities
shape_map <- c("Two"=16, "Four"=17, "Eight"=15)  # ● ▲ ■

theme_nat <- theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 7, colour = "black"),
    axis.line         = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks.length = unit(1.6, "mm"),
    legend.position   = "right",
    legend.title      = element_text(size = 8, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(3.5, "mm"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.title        = element_blank()
  )

choose_ok_stresses <- function(D, xcol) {
  D |>
    filter(is.finite(.data[[xcol]]), .data[[xcol]] > 0) |>
    group_by(Stress) |>
    summarise(span = max(.data[[xcol]], na.rm = TRUE) / min(.data[[xcol]], na.rm = TRUE), .groups = "drop") |>
    filter(is.finite(span), span >= 2) |>
    pull(Stress) |> as.character()
}

# ---------- UNWEIGHTED ----------
Du <- df %>% rename(x = Trait_var) %>% filter(is.finite(x), x > 0)
ok_u <- choose_ok_stresses(Du, "x")
dec_breaks_u <- 10^(seq(floor(log10(min(Du$x))), ceiling(log10(max(Du$x))), by = 1))

p_unw_nat <- ggplot(Du, aes(x = x, y = log10_Yield, colour = Stress, shape = Richness_f)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_smooth(
    data = subset(Du, Stress %in% ok_u),
    mapping = aes(x = x, y = log10_Yield, colour = Stress, group = Stress),
    inherit.aes = FALSE,
    method = "lm", se = FALSE, linewidth = 0.7
  ) +
  scale_x_log10(breaks = dec_breaks_u, labels = math_format(10^.x)) +
  annotation_logticks(sides = "b") +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    x = "Variance in response traits (growth rate)",
    y = expression(log[10]*"(Community yield)")
  ) +
  theme_nat

ggsave("plots/Function_vs_TraitVariance_UNW_Nature.tiff",
       p_unw_nat, width = 90, height = 70, units = "mm", dpi = 600, compression = "lzw")
ggsave("plots/Function_vs_TraitVariance_UNW_Nature.png",
       p_unw_nat, width = 90, height = 70, units = "mm", dpi = 600)

# ---------- WEIGHTED ----------
Dw <- df %>% rename(x = Trait_var_w) %>% filter(is.finite(x), x > 0)
ok_w <- choose_ok_stresses(Dw, "x")
dec_breaks_w <- 10^(seq(floor(log10(min(Dw$x))), ceiling(log10(max(Dw$x))), by = 1))

p_w_nat <- ggplot(Dw, aes(x = x, y = log10_Yield, colour = Stress, shape = Richness_f)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_smooth(
    data = subset(Dw, Stress %in% ok_w),
    mapping = aes(x = x, y = log10_Yield, colour = Stress, group = Stress),
    inherit.aes = FALSE,
    method = "lm", se = FALSE, linewidth = 0.7
  ) +
  scale_x_log10(breaks = dec_breaks_w, labels = math_format(10^.x)) +
  annotation_logticks(sides = "b") +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = shape_map, name = "Richness") +
  labs(
    x = "Abundance-weighted variance in response traits (growth rate)",
    y = expression(log[10]*"(Community yield)")
  ) +
  theme_nat

ggsave("plots/Function_vs_TraitVariance_WEIGHTED_Nature.tiff",
       p_w_nat, width = 90, height = 70, units = "mm", dpi = 600, compression = "lzw")
ggsave("plots/Function_vs_TraitVariance_WEIGHTED_Nature.png",
       p_w_nat, width = 90, height = 70, units = "mm", dpi = 600)

# ---------- Side-by-side panel ----------
ggpubr::ggarrange(p_unw_nat, p_w_nat, ncol = 2, labels = c("a","b")) |>
  ggpubr::ggexport(filename = "plots/Function_vs_TraitVariance_Nature_panels.png",
                   width = 1900, height = 800, res = 300)

message("Saved:\n- plots/Function_vs_TraitVariance_UNW_Nature.(tiff|png)\n- plots/Function_vs_TraitVariance_WEIGHTED_Nature.(tiff|png)\n- plots/Function_vs_TraitVariance_Nature_panels.png")



# ================================================================
# Extract stressor combinations from each CSV and combine to ONE CSV
# - Keeps zero estimates (no epsilon)
# - Normalizes levels to: Temp {20, 38}, Salinity {0, 20}, pH {7.2, 5.5}
# - Derives unified Stress label: Control, Temp, pH, Sal, pHTemp, SalTemp, pHSal, pHSalTemp
# - Outputs:
#     tables/Growth_all_sources_stressors_long.csv       (full rows)
#     tables/Growth_all_sources_stressors_summary.csv    (n, mean, sd per Id×Stress×Source)
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# -------------------- IO --------------------
IN_TEMP <- "LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv"
IN_pH   <- "LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv"
IN_SAL  <- "LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv"

OUT_DIR <- "tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------- Helpers --------------------
trimf    <- function(x) trimws(as.character(x))
safe_num <- function(x) suppressWarnings(as.numeric(x))
has_col  <- function(d, nm) nm %in% names(d)

sanitize_names <- function(df) {
  nm <- names(df)
  nm <- gsub("\\p{Zs}+", " ", nm, perl = TRUE)  # collapse weird spaces
  nm <- trimws(nm)
  for (i in seq_along(nm)) if (is.na(nm[i]) || nm[i] == "") nm[i] <- paste0("V", i)
  names(df) <- make.unique(nm, sep = "_")
  df
}

load_csv <- function(path){
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) |>
    sanitize_names()
}

# Standardise core columns (taxon & growth) without dropping zeros
standardise_core_cols <- function(d){
  id_col <- intersect(c("Id","OTU","Taxon","Strain","taxon","id"), names(d))[1]
  g_col  <- intersect(c("estimate","g","growth","rate","r","Growth","Rate"), names(d))[1]
  if (is.na(id_col)) stop("Could not find a taxon column among: Id/OTU/Taxon/Strain.")
  if (is.na(g_col))  stop("Could not find a growth column among: estimate/g/rate/r.")
  d %>%
    mutate(
      Id = trimf(.data[[id_col]]),
      g  = safe_num(.data[[g_col]])
    )
}

# Snap levels to discrete sets (others -> NA so we don't mislabel)
norm_temp <- function(x){
  x <- safe_num(x); x[!(x %in% c(20, 38))] <- NA_real_; x
}
norm_sal <- function(x){
  x <- safe_num(x); x[!(x %in% c(0, 20))]  <- NA_real_; x
}
norm_pH <- function(x){
  x <- safe_num(x)
  # keep 5.5 as 5.5, and treat 7.2 (or 7.0–7.4 if present) as 7.2
  case_when(
    is.na(x)            ~ NA_real_,
    abs(x - 5.5) <= 0.2 ~ 5.5,
    abs(x - 7.2) <= 0.3 ~ 7.2,
    TRUE                ~ NA_real_
  )
}

# Derive unified Stress label from normalized levels
derive_stress <- function(temp, sal, ph){
  temp <- norm_temp(temp); sal <- norm_sal(sal); ph <- norm_pH(ph)
  ifelse(is.na(temp) | is.na(sal) | is.na(ph), NA_character_,
         {
           t_on  <- ifelse(temp == 38, 1L, 0L)  # 38 vs 20
           s_on  <- ifelse(sal  == 20, 1L, 0L)  # 20 vs 0
           p_on  <- ifelse(ph   == 5.5, 1L, 0L) # 5.5 vs 7.2
           key <- paste(p_on, s_on, t_on, sep = "")
           dplyr::recode(key,
                         "000" = "Control",
                         "001" = "Temp",
                         "010" = "Sal",
                         "011" = "SalTemp",
                         "100" = "pH",
                         "101" = "pHTemp",
                         "110" = "pHSal",
                         "111" = "pHSalTemp",
                         .default = NA_character_
           )
         })
}

# Clean one file -> long tidy with Source + normalized levels + Stress
clean_one <- function(path, source_tag){
  raw <- load_csv(path) |> standardise_core_cols()
  raw %>%
    mutate(
      Source         = source_tag,
      TempLevel_raw  = if (has_col(., "TempLevel"))     safe_num(TempLevel)     else NA_real_,
      Salinity_raw   = if (has_col(., "SalinityLevel")) safe_num(SalinityLevel) else NA_real_,
      pH_raw         = if (has_col(., "pH_Level"))      safe_num(pH_Level)      else NA_real_,
      TempLevel      = norm_temp(TempLevel_raw),
      SalinityLevel  = norm_sal(Salinity_raw),
      pH_Level       = norm_pH(pH_raw),
      Stress         = derive_stress(TempLevel, SalinityLevel, pH_Level)
    ) %>%
    # keep rows that map cleanly to a discrete combination
    filter(!is.na(Stress)) %>%
    # keep core visibility; retain Rep if present, otherwise NA
    mutate(Rep = if (has_col(., "Rep")) as.integer(Rep) else as.integer(NA)) %>%
    select(Source, Id, g, Rep, TempLevel, SalinityLevel, pH_Level, Stress)
}

# -------------------- Run for all three inputs --------------------
d_temp <- clean_one(IN_TEMP, "Temp")
d_pH   <- clean_one(IN_pH,   "pH")
d_sal  <- clean_one(IN_SAL,  "Sal")

all_long <- bind_rows(d_temp, d_pH, d_sal) %>%
  mutate(
    Id     = trimf(Id),
    Stress = factor(Stress, levels = c("Control","Temp","pH","Sal",
                                       "pHTemp","SalTemp","pHSal","pHSalTemp"))
  ) %>%
  # Keep non-negative g; zeros are real and kept. Drop NA g rows only.
  filter(!is.na(g), g >= 0)

# Write ONE combined CSV
write.csv(all_long,
          file.path(OUT_DIR, "Growth_all_sources_stressors_long.csv"),
          row.names = FALSE)

# Also write a compact summary table (per Source × Id × Stress)
summary_tbl <- all_long %>%
  group_by(Source, Id, Stress) %>%
  summarise(
    n      = n(),
    mean_g = mean(g),
    sd_g   = sd(g),
    .groups = "drop"
  ) %>%
  arrange(Source, Id, Stress)

write.csv(summary_tbl,
          file.path(OUT_DIR, "Growth_all_sources_stressors_summary.csv"),
          row.names = FALSE)

message("Wrote:\n - ", file.path(OUT_DIR, "Growth_all_sources_stressors_long.csv"),
        "\n - ", file.path(OUT_DIR, "Growth_all_sources_stressors_summary.csv"))


suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(rstatix)
  library(glue)
  library(scales)
})

# ----------------------
# Helpers
# ----------------------
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  ggplot2::scale_x_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# ----------------------
# User parameters
# ----------------------
PATH_CSV      <- "Growth_all_sources_stressors_long.csv"
MIN_REPS      <- 2
P_ADJ_METHOD  <- "BH"
OUT_DIR       <- "FOCAL_vs_OTHERS_8FACET_HEATMAP"

# ----------------------
# Setup & data
# ----------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(PATH_CSV, show_col_types = FALSE) %>%
  clean_names() %>%
  filter(!is.na(id), !is.na(g), !is.na(stress)) %>%
  mutate(
    id = as.character(id),
    stress = as.factor(stress),
    group = substr(id, 1, 1) # e.g., "I" or "D"
  )

all_taxa <- sort(unique(df$id))
if (length(all_taxa) < 2) stop("Need at least 2 taxa for pairwise comparisons.")

# ----------------------
# Pairwise stats (focal vs other) per Stress
# ----------------------
stress_levels <- levels(droplevels(df$stress))
res_list <- vector("list", 0L)

for (fid in all_taxa) {
  for (st in stress_levels) {
    sub <- df %>% filter(stress == st)
    if (n_distinct(sub$id) < 2) next
    others <- setdiff(unique(sub$id), fid)
    if (length(others) == 0) next
    
    for (oth in others) {
      A <- sub %>% filter(id == fid)
      B <- sub %>% filter(id == oth)
      n_a <- nrow(A); n_b <- nrow(B)
      
      if ((n_a < MIN_REPS) | (n_b < MIN_REPS)) {
        res_list[[length(res_list) + 1]] <- tibble(
          focal = fid, other = oth, stress = st,
          group_focal = substr(fid, 1, 1), group_other = substr(oth, 1, 1),
          n_focal = n_a, n_other = n_b,
          mean_focal = mean(A$g, na.rm = TRUE), mean_other = mean(B$g, na.rm = TRUE),
          diff_mean = mean(A$g, na.rm = TRUE) - mean(B$g, na.rm = TRUE),
          effect_size = NA_real_, effect_ci_low = NA_real_, effect_ci_high = NA_real_,
          test = NA_character_, p_value = NA_real_, sufficient_reps = FALSE
        )
        next
      }
      
      pair_dat <- bind_rows(
        A %>% mutate(comp_group = "FOCAL"),
        B %>% mutate(comp_group = "OTHER")
      )
      
      t_res <- tryCatch(
        t_test(g ~ comp_group, data = pair_dat, var.equal = FALSE),
        error = function(e) NULL
      )
      es <- tryCatch(
        cohens_d(g ~ comp_group, data = pair_dat, pooled_sd = FALSE, hedges_correction = TRUE),
        error = function(e) NULL
      )
      
      pval <- if (!is.null(t_res) && nrow(t_res) == 1) t_res$p[1] else NA_real_
      eff  <- if (!is.null(es) && nrow(es) == 1) es$effsize[1] else NA_real_
      ci_l <- if (!is.null(es) && nrow(es) == 1) es$conf.low[1] else NA_real_
      ci_h <- if (!is.null(es) && nrow(es) == 1) es$conf.high[1] else NA_real_
      
      res_list[[length(res_list) + 1]] <- tibble(
        focal = fid, other = oth, stress = st,
        group_focal = substr(fid, 1, 1), group_other = substr(oth, 1, 1),
        n_focal = n_a, n_other = n_b,
        mean_focal = mean(A$g, na.rm = TRUE), mean_other = mean(B$g, na.rm = TRUE),
        diff_mean = mean(A$g, na.rm = TRUE) - mean(B$g, na.rm = TRUE),
        effect_size = eff, effect_ci_low = ci_l, effect_ci_high = ci_h,
        test = "Welch t", p_value = pval, sufficient_reps = TRUE
      )
    }
  }
}

pairwise <- bind_rows(res_list)
if (nrow(pairwise) == 0) stop("No pairwise results.")

pairwise <- pairwise %>%
  mutate(p_adj = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = P_ADJ_METHOD)))

write_csv(pairwise, file.path(OUT_DIR, "pairwise_results.csv"))

# ----------------------
# 8-facet HEATMAP (facet by Stress only) with significance asterisks
# ----------------------
# nice Stress order (customize to your set)
nice_stress_order <- c("Control", "pH", "Sal", "Temp", "pHSal", "pHTemp", "SalTemp", "pHSalTemp")

# global row/col order to keep matrices aligned across facets
row_order <- sort(unique(pairwise$focal))
col_order <- sort(unique(pairwise$other))

hm <- pairwise %>%
  mutate(
    stress = factor(as.character(stress),
                    levels = unique(c(nice_stress_order,
                                      setdiff(sort(unique(stress)), nice_stress_order)))),
    focal = factor(focal, levels = row_order),
    other = factor(other, levels = col_order),
    sig_lab = case_when(
      is.na(p_adj) | !sufficient_reps ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(focal != other)  # remove diagonal

max_abs <- max(abs(hm$diff_mean), na.rm = TRUE)

p_heat8 <- ggplot(hm, aes(x = other, y = focal, fill = diff_mean)) +
  geom_tile() +
  geom_text(aes(label = sig_lab), fontface = "bold", size = 3, na.rm = TRUE) +
  scale_fill_distiller(
    palette = "RdBu", direction = 1,
    limits = c(-max_abs, max_abs), oob = squish,
    na.value = "grey92",
    name = "Δ mean g\n(focal − other)"
  ) +
  facet_wrap(~ stress, ncol = 4, scales = "free") +
  labs(
    title = "Relative growth (Δ mean g): focal vs other",
    subtitle = "8-facet heatmap by Stress. Asterisks denote BH-adjusted significance (*** <0.001, ** <0.01, * <0.05).",
    x = "Other taxon", y = "Focal taxon"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "HEATMAP_8facet_asterisks.png"), p_heat8, width = 14, height = 10, dpi = 300)

# ----------------------
# Holistic DOT plot — z-scored within Stress, short labels only
# ----------------------

# id -> short label (no IDs)
taxon_map <- tibble::tribble(
  ~id,  ~id_lab,
  "D11","Ped.",
  "D14","Mic.",
  "D17","Cur.",
  "I2", "Yer.",
  "I8", "Pse.F.",
  "I9", "Chr.",
  "I11","Pse.A.",
  "I15","Ser.",
  "I18","Aer.",
  "I20","Chry.",
  "I22","Erw.",
  "I23","Jan."
)

df_lab <- df %>%
  left_join(taxon_map, by = "id") %>%
  mutate(id_lab = if_else(is.na(id_lab), id, id_lab))  # fallback if any missing

# z-score g within each Stress
z_df <- df_lab %>%
  group_by(stress) %>%
  mutate(g_z = as.numeric(scale(g))) %>%
  ungroup()

# mean z per stress × taxon (for ordering)
mean_by_stress_taxon <- z_df %>%
  group_by(stress, id_lab) %>%
  summarise(mu = mean(g_z, na.rm = TRUE), .groups = "drop")

plot_df <- z_df %>%
  left_join(mean_by_stress_taxon, by = c("stress","id_lab"))

p_dot_short <- ggplot(plot_df, aes(x = reorder_within(id_lab, mu, stress), y = g_z)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 1.8, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.4, color = "black", alpha = 0.95) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.28, color = "black", alpha = 0.95) +
  facet_wrap(~ stress, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  labs(
    title = "Standardized growth (z within Stress) by taxon",
    x = "Taxon (ordered within facet)", y = "z-scored g"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(face = "italic")   # ← italicize genus names
  )


ggsave(file.path(OUT_DIR, "DOTPLOT_z_by_stress_shortlabels.png"),
       p_dot_short, width = 12, height = 7, dpi = 300)


# ----------------------
# Compact DOT plot — z-scored within Stress
#  • short italic labels
#  • faint replicates
#  • colored means (diverging by mean z)
#  • tight facet spacing & minimal grid
# ----------------------

# id -> short label (no IDs)
taxon_map <- tibble::tribble(
  ~id,  ~id_lab,
  "D11","Ped.", "D14","Mic.", "D17","Cur.", "I2","Yer.",
  "I8","Pse.F.", "I9","Chr.", "I11","Pse.A.", "I15","Ser.",
  "I18","Aer.", "I20","Chry.", "I22","Erw.", "I23","Jan."
)

df_lab <- df %>%
  left_join(taxon_map, by = "id") %>%
  mutate(id_lab = if_else(is.na(id_lab), id, id_lab))  # fallback if any missing

# z-score within each Stress
z_df <- df_lab %>%
  group_by(stress) %>%
  mutate(g_z = as.numeric(scale(g))) %>%
  ungroup()

# summary for ordering + coloring
sumz <- z_df %>%
  group_by(stress, id_lab) %>%
  summarise(
    mu = mean(g_z, na.rm = TRUE),
    se = sd(g_z, na.rm = TRUE) / sqrt(sum(is.finite(g_z))),
    .groups = "drop"
  )

plot_df <- z_df %>% left_join(sumz, by = c("stress","id_lab"))

# symmetric limits per facet for color balance (optional: global)
lim_global <- max(abs(sumz$mu), na.rm = TRUE)

library(ggplot2)
p_dot_compact <- ggplot() +
  # faint replicates
  geom_point(
    data = plot_df,
    aes(x = reorder_within(id_lab, mu, stress), y = g_z),
    size = 1.1, alpha = 0.35, color = "grey30",
    position = position_jitter(width = 0.08, height = 0)
  ) +
  # mean ± SE (thin)
  geom_errorbar(
    data = sumz,
    aes(x = reorder_within(id_lab, mu, stress),
        ymin = mu - se, ymax = mu + se),
    width = 0.18, linewidth = 0.3, color = "black"
  ) +
  # colored mean point
  geom_point(
    data = sumz,
    aes(x = reorder_within(id_lab, mu, stress), y = mu, fill = mu),
    shape = 21, size = 2.8, stroke = 0.3, color = "black"
  ) +
  # helpful reference lines
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey40") +
  facet_wrap(~ stress, scales = "free", ncol = 4) +
  coord_flip(clip = "off") +
  scale_x_reordered() +
  # diverging fill for means (consistent across facets)
  scale_fill_distiller(
    palette = "RdBu", direction = 1,
    limits = c(-lim_global, lim_global), oob = scales::squish,
    name = "z-score",   # label for legend
    guide = guide_colorbar(
      title.position = "right",  # put title next to bar
      title.hjust = 0.5
    )
  ) +
  labs(
    x = NULL, y = "z-scored g"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88"),
    axis.text.y = element_text(face = "italic", size = 8, margin = margin(r = 2)),
    axis.text.x = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "top",    # legend at top
    legend.direction = "horizontal",
    panel.spacing.x = unit(4, "pt"),
    panel.spacing.y = unit(4, "pt"),
    plot.margin = margin(6, 8, 6, 6)
  )

ggsave(file.path(OUT_DIR, "DOTPLOT_z_by_stress_compact.png"),
       p_dot_compact, width = 11, height = 6.2, dpi = 300)

# export replicate-level z-scores
readr::write_csv(z_df, file.path(OUT_DIR, "Zscores_replicates.csv"))

# export summary means & SE
readr::write_csv(sumz, file.path(OUT_DIR, "Zscores_summary.csv"))




# ================================================================
# Bayesian softmax composition (Dirichlet or Dirichlet–multinomial)
# + Uncertainty-aware biomass (OD) via MI from posterior AWM draws
#
# Fully refined & reordered version
# - Allows positive/negative slopes for kappa (no exp)
# - Optional taxon-bias term (delta) gated via data flag
# - Clean DM model (built-in lpmf) + log_lik for LOO
# - Reproducible draw subsampling
# - Optional nested-CV model selection for OD mapping
# - Safer IO (file names), tidy helpers, PPC stubs
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(grid)
  library(rstan)
  library(loo)
  library(philentropy)   # JSD
  library(ragg)
})

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(), stringsAsFactors = FALSE)
set.seed(123)

# ---------------------- USER OPTIONS ---------------------------------
USE_RICHNESS_KAPPA   <- FALSE   # TRUE => kappa by Stress x Richness
USE_TAXON_BIASES     <- TRUE    # include delta_i with shrinkage prior
LIKELIHOOD            <- "dirichlet"  # "dirichlet" or "dm" (auto-switch later)

# OD model options
USE_NESTED_CV_FOR_OD <- FALSE   # set TRUE to select OD formula by nested CV
K_OUTER <- 5
K_INNER <- 3

# Stan sampling options
CHAINS   <- 4
ITER     <- 3000
WARMUP   <- 1500
ADAPT    <- list(adapt_delta = 0.99, max_treedepth = 15)

# Posterior summarisation
NDRAWS_COMPOSITION <- 400   # draws for composition summaries (p_hat)
NDRAWS_MI          <- 400   # draws for multiple-imputation AWM -> OD

# ---------------------- HELPERS --------------------------------------
trimf   <- function(x) trimws(as.character(x))
softmax <- function(x){ x <- x - max(x); ex <- exp(x); ex/sum(ex) }
recode_stress <- function(x){
  dplyr::case_when(
    x == "pH_Sal"      ~ "pHSal",
    x == "pH_Sal_Temp" ~ "pHSalTemp",
    x == "pH_Temp"     ~ "pHTemp",
    x == "Sal_Temp"    ~ "SalTemp",
    TRUE ~ x
  )
}

# ---------------------- DESIGN (planned members) ---------------------
design <- tribble(
  ~Com_Id, ~Id,
  "Com1","D14", "Com1","I15",
  "Com2","D17", "Com2","I9",
  "Com3","I8",  "Com3","I11",
  "Com4","I2",  "Com4","I18",
  "Com5","I20", "Com5","I23",
  "Com6","I22", "Com6","D11",
  # 4-taxon
  "Com7","I20","Com7","I22","Com7","I23","Com7","D14",
  "Com8","I2","Com8","I9","Com8","I15","Com8","D17",
  "Com9","I8","Com9","I11","Com9","I18","Com9","D14",
  "Com10","I2","Com10","I8","Com10","I15","Com10","I23",
  "Com11","I9","Com11","I11","Com11","I20","Com11","D11",
  "Com12","I18","Com12","I22","Com12","D11","Com12","D17",
  # 8-taxon
  "Com13","I2","Com13","I9","Com13","I11","Com13","I20",
  "Com13","I22","Com13","D11","Com13","D14","Com13","D17",
  "Com14","I2","Com14","I8","Com14","I9","Com14","I18",
  "Com14","I20","Com14","I22","Com14","D11","Com14","D17",
  "Com15","I9","Com15","I11","Com15","I15","Com15","I18",
  "Com15","I20","Com15","I23","Com15","D11","Com15","D14",
  "Com16","I11","Com16","I18","Com16","I20","Com16","I22",
  "Com16","I23","Com16","D11","Com16","D14","Com16","D17",
  "Com17","I2","Com17","I8","Com17","I9","Com17","I15",
  "Com17","I20","Com17","I23","Com17","D11","Com17","D17",
  "Com18","I2","Com18","I9","Com18","I15","Com18","I18",
  "Com18","I20","Com18","I22","Com18","D11","Com18","D17",
  "Com19","I2","Com19","I8","Com19","I9","Com19","I15",
  "Com19","I18","Com19","I20","Com19","I23","Com19","D14",
  "Com20","I2","Com20","I8","Com20","I9","Com20","I11",
  "Com20","I18","Com20","I20","Com20","I22","Com20","D14",
  "Com21","I2","Com21","I8","Com21","I11","Com21","I15",
  "Com21","I18","Com21","I20","Com21","I23","Com21","D14",
  "Com22","I8","Com22","I11","Com22","I15","Com22","I20",
  "Com22","I22","Com22","I23","Com22","D11","Com22","D17"
) %>% mutate(Com_Id = gsub("\\s+","", Com_Id), Id = trimf(Id))

# ---------------------- LOAD & PREP DATA -----------------------------
message("[1/7] Loading abundance data ...")

abund_long <- read.csv("All_Abundance_Data.csv") %>%  # fixed filename case
  transmute(
    SampleID  = trimf(SampleID),
    Com_Id    = gsub("\\s+", "", trimf(Com_Id)),
    Stress    = trimf(Stress),
    Diversity = trimf(Diversity),
    Id        = trimf(Id),
    Abundance = as.numeric(Abundance)
  ) %>%
  mutate(
    Diversity = dplyr::recode(Diversity,
                              "One"="1","one"="1","Two"="2","two"="2",
                              "Four"="4","four"="4","Eight"="8","eight"="8",
                              .default = Diversity),
    Diversity = as.numeric(Diversity),
    Stress    = recode_stress(Stress)
  ) %>%
  filter(Diversity %in% c(2,4,8)) %>%
  group_by(SampleID) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# Expand each sample to planned membership and re-normalize
abund_complete <- abund_long %>%
  distinct(SampleID, Com_Id, Stress, Diversity) %>%
  left_join(design, by = "Com_Id", relationship = "many-to-many") %>%
  left_join(abund_long, by = c("SampleID","Com_Id","Stress","Diversity","Id"),
            relationship = "many-to-many") %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  group_by(SampleID) %>% mutate(Abundance = Abundance/sum(Abundance)) %>% ungroup()

# ---------------------- GROWTH FILES & PREP --------------------------
message("[2/7] Loading growth rates ...")
load_growth_file <- function(path, mutate_steps){ mutate_steps(read.csv(path)) }

mut_temp <- function(d){ d %>%
    filter(TempLevel %in% c(20, 38)) %>%
    mutate(Stress = if_else(TempLevel==38, paste0(Stress, "Temp"), Stress)) %>%
    mutate(Stress = case_when(Stress=="ControlTemp"~"Temp",
                              Stress=="pH_SalTemp"~"pHSalTemp",
                              Stress=="pH_Sal"~"pHSal",
                              TRUE ~ Stress)) %>%
    select(-TempLevel, -Rep)
}
mut_pH <- function(d){ d %>%
    mutate(Stress = case_when((SalinityLevel==0 & TempLevel==20)~"Control",
                              (SalinityLevel==0 & TempLevel==38)~"Temp",
                              (SalinityLevel==20 & TempLevel==20)~"Sal",
                              (SalinityLevel==20 & TempLevel==38)~"SalTemp")) %>%
    select(-SalinityLevel, -TempLevel, -X, -X.1) %>%
    filter(pH_Level %in% c(5.5,7.2)) %>%
    mutate(Stress = if_else(pH_Level==5.5, paste0("pH", Stress), Stress)) %>%
    mutate(Stress = if_else(Stress=="pHControl","pH",Stress)) %>%
    select(-pH_Level, -Rep)
}
mut_sal <- function(d){ d %>%
    mutate(Stress = case_when((pH_Level==7.2 & TempLevel==20)~"Control",
                              (pH_Level==7.2 & TempLevel==38)~"Temp",
                              (pH_Level==5.5 & TempLevel==20)~"pH",
                              (pH_Level==5.5 & TempLevel==38)~"pHTemp")) %>%
    select(-pH_Level, -TempLevel, -X) %>%
    filter(SalinityLevel %in% c(0,20)) %>%
    mutate(Stress = if_else(SalinityLevel==20, paste0(Stress, "Sal"), Stress)) %>%
    mutate(Stress = case_when(Stress=="ControlSal"~"Sal",
                              Stress=="pHTempSal"~"pHSalTemp",
                              Stress=="TempSal"~"SalTemp",
                              TRUE ~ Stress)) %>%
    select(-SalinityLevel, -Rep)
}

all_growth_raw <- bind_rows(
  load_growth_file("LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv", mut_temp),
  load_growth_file("LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv", mut_pH),
  load_growth_file("LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv", mut_sal)
)

all_growth <- all_growth_raw %>%
  mutate(Stress = recode_stress(Stress)) %>%
  group_by(Id, Stress) %>%
  summarise(g = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  mutate(Id = trimf(Id)) %>%
  group_by(Stress) %>%
  mutate(g_z = as.numeric(scale(g))) %>%
  ungroup()

# ---------------------- LINK GROWTH & FILTER -------------------------
message("[3/7] Linking growth to abundance & filtering complete samples ...")

abund_growth <- abund_complete %>% left_join(all_growth, by = c("Id","Stress"))
complete_ids <- abund_growth %>% group_by(SampleID) %>%
  summarise(all_g = all(!is.na(g)), .groups = "drop") %>% filter(all_g) %>% pull(SampleID)

data_all <- abund_growth %>% filter(SampleID %in% complete_ids)

# ---------------------- CHOOSE LIKELIHOOD ----------------------------
if (LIKELIHOOD == "dm" && !("Reads" %in% names(data_all))) {
  warning("LIKELIHOOD='dm' requested but no 'Reads' found; using 'dirichlet'.")
  LIKELIHOOD <- "dirichlet"
}

# ---------------------- RAGGED STRUCTURES ----------------------------
message("[4/7] Building ragged sample structures ...")

ord <- data_all %>% arrange(SampleID, Id)
SAMPLES <- ord %>% group_by(SampleID) %>% summarise(
  Com_Id    = first(Com_Id),
  Stress    = first(Stress),
  Diversity = first(Diversity),
  len       = n(),
  .groups   = "drop"
)

N  <- nrow(SAMPLES)
len <- SAMPLES$len
start_idx <- c(1L, 1L + head(cumsum(len), -1L))
J  <- sum(len)

p_obs  <- ord$Abundance
g_z    <- ord$g_z
g_raw  <- ord$g

K <- ord %>% distinct(Id) %>% arrange(Id) %>% nrow()
taxa_levels <- ord %>% distinct(Id) %>% arrange(Id) %>% pull(Id)
map_tax <- setNames(seq_along(taxa_levels), taxa_levels)

tax_of      <- as.integer(map_tax[ord$Id])
stress_lvls <- sort(unique(ord$Stress))
S <- length(stress_lvls)
map_stress <- setNames(seq_along(stress_lvls), stress_lvls)
stress_id   <- as.integer(map_stress[SAMPLES$Stress])

rich_lvls <- sort(unique(SAMPLES$Diversity))
Rdim <- if (USE_RICHNESS_KAPPA) length(rich_lvls) else 1L
rich_id <- if (USE_RICHNESS_KAPPA) as.integer(factor(SAMPLES$Diversity, levels = rich_lvls)) else rep(1L, N)

if (LIKELIHOOD == "dm") {
  if (!("Reads" %in% names(ord))) stop("DM selected but 'Reads' column not in 'ord'.")
  y <- as.integer(ord$Reads)
  if (any(is.na(y))) stop("NA in Reads; supply true counts or switch to Dirichlet.")
}

# sanity: each segment sums to 1
seg_chk <- purrr::map_dfr(seq_len(N), function(n){
  a <- start_idx[n]; L <- len[n]
  tibble(n=n, sum_seg = sum(p_obs[a:(a+L-1)]))
})
stopifnot(all(abs(seg_chk$sum_seg - 1) < 1e-8))

# ---------------------- STAN MODELS ----------------------------------
message("[5/7] Writing Stan model(s) ...")

stan_dirichlet <- '
functions {
  vector segment_vector(vector x, int start, int L) {
    vector[L] y;
    for (m in 1:L) y[m] = x[start + m - 1];
    return y;
  }
}
data {
  int<lower=1> N;                 // samples
  int<lower=1> K;                 // total taxa
  int<lower=1> J;                 // total elements in ragged arrays
  int<lower=1> S;                 // stresses
  int<lower=1> start_idx[N];      // start pos for sample n
  int<lower=1> len[N];            // length for sample n
  vector[J] p_obs;                // observed proportions (flattened)
  vector[J] g_z;                  // standardized growth for logits
  int<lower=1, upper=K> tax_of[J];// taxon index per element
  int<lower=1, upper=S> stress_id[N];
  int<lower=1> Rdim;              // richness dimension actually used
  int<lower=1, upper=Rdim> rich_id[N];
  int<lower=0, upper=1> use_taxon_biases; // gate delta
}
parameters {
  // kappa effects (non-centered)
  real mu_kappa;
  real<lower=0> sigma_kappa;
  matrix[S, Rdim] kappa_raw;

  // taxon biases
  vector[K-1] delta_raw;
  real<lower=0> sigma_delta;

  // log-phi hierarchy
  real mu_phi;
  real<lower=0> sigma_phi;
  vector[S] log_phi_raw;
}
transformed parameters {
  matrix[S, Rdim] kappa; // can be +/-
  vector[K] delta;
  vector[S] log_phi;
  vector[S] phi;

  kappa = mu_kappa + sigma_kappa * kappa_raw;

  delta[1:(K-1)] = use_taxon_biases * sigma_delta * delta_raw;
  delta[K] = 0;

  log_phi = mu_phi + sigma_phi * log_phi_raw;
  for (s in 1:S) phi[s] = exp(log_phi[s]);
}
model {
  // priors
  mu_kappa     ~ normal(0, 1);
  sigma_kappa  ~ normal(0, 0.5);
  to_vector(kappa_raw) ~ normal(0, 1);

  sigma_delta  ~ normal(0, 1);
  delta_raw    ~ normal(0, 1);

  mu_phi       ~ normal(log(50), 1);
  sigma_phi    ~ normal(0, 1);
  log_phi_raw  ~ normal(0, 1);

  // likelihood: Dirichlet on observed proportions
  for (n in 1:N) {
    int a = start_idx[n];
    int L = len[n];
    vector[L] eta;
    vector[L] p_n;
    vector[L] pseg = segment_vector(p_obs, a, L);

    // tiny jitter + renorm for any exact zeros
    pseg = pseg + 1e-12;
    pseg = pseg / sum(pseg);

    for (m in 1:L) {
      int j = a + m - 1;
      int t = tax_of[j];
      eta[m] = kappa[stress_id[n], rich_id[n]] * g_z[j] + delta[t];
    }
    p_n = softmax(eta);
    target += dirichlet_lpdf(pseg | phi[stress_id[n]] * p_n);
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

stan_dm <- '
functions {
  int[] segment_int_array(int[] x, int start, int L) {
    int y[L];
    for (m in 1:L) y[m] = x[start + m - 1];
    return y;
  }
}
data {
  int<lower=1> N; int<lower=1> K; int<lower=1> J; int<lower=1> S;
  int<lower=1> start_idx[N]; int<lower=1> len[N];
  int<lower=0> y[J]; vector[J] g_z;
  int<lower=1, upper=K> tax_of[J]; int<lower=1, upper=S> stress_id[N];
  int<lower=1> Rdim; int<lower=1, upper=Rdim> rich_id[N];
  int<lower=0, upper=1> use_taxon_biases;
}
parameters {
  real mu_kappa; real<lower=0> sigma_kappa; matrix[S, Rdim] kappa_raw;
  vector[K-1] delta_raw; real<lower=0> sigma_delta;
  real mu_phi; real<lower=0> sigma_phi; vector[S] log_phi_raw;
}
transformed parameters {
  matrix[S, Rdim] kappa;
  vector[K] delta; vector[S] log_phi; vector[S] phi;
  kappa = mu_kappa + sigma_kappa * kappa_raw;
  delta[1:(K-1)] = use_taxon_biases * sigma_delta * delta_raw; delta[K] = 0;
  log_phi = mu_phi + sigma_phi * log_phi_raw; for (s in 1:S) phi[s] = exp(log_phi[s]);
}
model {
  mu_kappa ~ normal(0,1); sigma_kappa ~ normal(0,0.5); to_vector(kappa_raw) ~ normal(0,1);
  sigma_delta ~ normal(0,1); delta_raw ~ normal(0,1);
  mu_phi ~ normal(log(50),1); sigma_phi ~ normal(0,1); log_phi_raw ~ normal(0,1);
  for (n in 1:N) {
    int a = start_idx[n]; int L = len[n]; vector[L] eta; int seg[L] = segment_int_array(y, a, L);
    for (m in 1:L) {
      int j = a + m - 1; int t = tax_of[j];
      eta[m] = kappa[stress_id[n], rich_id[n]] * g_z[j] + delta[t];
    }
    target += dirichlet_multinomial_lpmf(seg | phi[stress_id[n]] * softmax(eta));
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    int a = start_idx[n]; int L = len[n]; vector[L] eta; int seg[L] = segment_int_array(y, a, L);
    for (m in 1:L) {
      int j = a + m - 1; int t = tax_of[j];
      eta[m] = kappa[stress_id[n], rich_id[n]] * g_z[j] + delta[t];
    }
    log_lik[n] = dirichlet_multinomial_lpmf(seg | phi[stress_id[n]] * softmax(eta));
  }
}
'

# Write model to disk and compile
if (LIKELIHOOD == "dirichlet") {
  writeLines(stan_dirichlet, "softmax_dirichlet.stan")
  sm <- stan_model("softmax_dirichlet.stan")
} else {
  writeLines(stan_dm, "softmax_dm.stan")
  sm <- stan_model("softmax_dm.stan")
}

# ---------------------- STAN DATA LIST -------------------------------
stan_data <- list(N=N, K=K, J=J, S=S,
                  start_idx=start_idx, len=len,
                  g_z=g_z, tax_of=tax_of, stress_id=stress_id,
                  Rdim=Rdim, rich_id=rich_id,
                  use_taxon_biases = as.integer(USE_TAXON_BIASES))

if (LIKELIHOOD == "dirichlet") stan_data$p_obs <- as.vector(p_obs)
if (LIKELIHOOD == "dm")        stan_data$y     <- if (exists("y")) y else stop("DM needs y")

# ---------------------- SAMPLE --------------------------------------
message("[6/7] Sampling ...")
fit <- sampling(sm, data = stan_data, chains = CHAINS, iter = ITER, warmup = WARMUP,
                control = ADAPT, seed = 123)

# Diagnostics
rstan::check_hmc_diagnostics(fit)
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
sum_div <- sum(sapply(sp, function(x) sum(x[,"divergent__"])))
cat("Total divergences:", sum_div, "\n")
ss <- summary(fit)$summary
cat("Max Rhat:", max(ss[,"Rhat"], na.rm=TRUE), "\n")
cat("Min ESS:", min(ss[,"n_eff"], na.rm=TRUE), "\n")

# Save bundle for reproducibility
saveRDS(list(fit=fit, stan_data=stan_data, SAMPLES=SAMPLES, ord=ord,
             meta=list(LIKELIHOOD=LIKELIHOOD, USE_TAXON_BIASES=USE_TAXON_BIASES,
                       USE_RICHNESS_KAPPA=USE_RICHNESS_KAPPA)),
        file = ifelse(LIKELIHOOD=="dirichlet","fit_dirichlet_bundle.rds","fit_dm_bundle.rds"))

# Print key parameters
print(fit, pars = c("kappa","sigma_kappa","mu_phi","sigma_phi"), probs = c(0.05,0.5,0.95))

# ---------------------- LOO -----------------------------------------
log_lik <- rstan::extract(fit, pars = "log_lik")$log_lik
if (!is.null(log_lik)) {
  loo_res <- loo::loo(log_lik)
  print(loo_res)
}

# ================== POSTERIOR COMPOSITIONS, AWM, METRICS & PLOTS ==================
message("[Post] Posterior compositions, AWM & evaluation ...")

stopifnot(exists("fit"), exists("SAMPLES"), exists("len"), exists("J"), exists("N"))
stopifnot(exists("tax_of"), exists("stress_id"))

# Extract posterior draws
draws <- rstan::extract(fit)
nd <- if (!is.null(draws$kappa)) dim(draws$kappa)[1] else length(draws$mu_kappa)

# Reproducible subsampling of draws
set.seed(123)
sel    <- if (nd > NDRAWS_COMPOSITION) sample.int(nd, NDRAWS_COMPOSITION) else seq_len(nd)
set.seed(123)
sel_mi <- if (nd > NDRAWS_MI)          sample.int(nd, NDRAWS_MI)          else seq_len(nd)

# Build segment indices once
idx_split <- split(seq_len(J), rep.int(seq_len(N), times = len))

# Growth vectors aligned with flattened order
g_z_vec   <- as.numeric(g_z)
# For AWM we want RAW growth (across stresses comparable); use ord$g
g_raw_vec <- as.numeric(g_raw)

# Helper to get kappa for a given draw/sample
get_kappa <- function(d, s, r = 1L) {
  if (length(dim(draws$kappa)) == 2L) {
    draws$kappa[d, s]
  } else if (length(dim(draws$kappa)) == 3L) {
    draws$kappa[d, s, r]
  } else stop("Unexpected 'kappa' dims")
}

# Posterior summaries per sample → per taxon; and AWM per sample
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
    d <- sel[ii]
    kap <- get_kappa(d, s, r_idx)
    logits <- kap * gsub_z + draws$delta[d, taxa_ix]
    p <- softmax(logits)
    P[ii, ] <- p
    A[ii]   <- sum(p * gsub_r)        # AWM with *raw* growth
  }
  
  p_mean <- colMeans(P)
  p_q    <- apply(P, 2, stats::quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
  
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

saveRDS(pred_comp, "bayes_pred_comp.rds")
saveRDS(AWM_post,  "bayes_awms.rds")

# Join observed compositions
top_obs <- ord %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(p_obs = Abundance) %>%
  dplyr::ungroup() %>%
  dplyr::select(SampleID, Id, p_obs)

pred_comp <- pred_comp %>% dplyr::left_join(top_obs, by = c("SampleID","Id"))

# ===================== METRICS =====================
comp_eval <- pred_comp %>%
  dplyr::filter(!is.na(p_obs)) %>%
  dplyr::summarise(
    RMSE = sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE)),
    R2   = 1 - sum((p_obs - p_hat)^2, na.rm = TRUE) /
      sum((p_obs - mean(p_obs, na.rm = TRUE))^2, na.rm = TRUE)
  )
print(comp_eval)

# --------------------- Composition PPC (quick) -----------------------
# For a random subset of samples, simulate p~Dir(phi*p_hat) and compare JSD
set.seed(123)
ppc_samples <- sample(unique(pred_comp$SampleID), size = min(12, length(unique(pred_comp$SampleID))))
ppc_df <- pred_comp %>% filter(SampleID %in% ppc_samples)
ppc_stats <- ppc_df %>% group_by(SampleID, Stress, Diversity) %>% summarise(
  JS_obs = philentropy::JSD(rbind(p_obs, p_hat)), .groups = "drop"
)
write.csv(ppc_stats, "PPC_composition_JS.csv", row.names = FALSE)

# ===================== PLOTS (with 90% CI) ===========================
comp_plot_df <- pred_comp %>% filter(!is.na(p_obs)) %>% mutate(Diversity = factor(Diversity))

comp_per_sample <- comp_plot_df %>%
  group_by(SampleID, Stress, Diversity) %>%
  summarise(
    RMSE = sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE)),
    JS   = philentropy::JSD(rbind(p_obs, p_hat)),
    .groups = "drop"
  )

lab_comp <- comp_per_sample %>%
  group_by(Stress) %>%
  summarise(
    RMSE = mean(RMSE, na.rm = TRUE),
    JS   = mean(JS,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("RMSE = ", sprintf("%.3f", RMSE),
                   "\nJS = ",   sprintf("%.3f", JS)),
    x = Inf, y = Inf
  )

p_comp <- ggplot(comp_plot_df, aes(x = p_hat, y = p_obs, colour = Diversity)) +
  geom_abline(linetype = 2) +
  geom_errorbarh(aes(xmin = p_lwr, xmax = p_upr), alpha = 0.30, height = 0, linewidth = 0.3) +
  geom_point(alpha = 0.85, size = 1.8) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  facet_wrap(~ Stress, nrow = 2) +
  scale_colour_manual(values = c("2"="#56B4E9","4"="#F0E442","8"="grey40"), name = "Richness") +
  labs(x = "Predicted relative abundance (posterior mean ± 90% CI)", y = "Observed relative abundance") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black", fill=NA, linewidth=0.4), panel.spacing = grid::unit(1.5, "lines")) +
  geom_text(data = lab_comp, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3)

ggsave("FIG_BAYES_comp_per_stressor.tiff", p_comp, width = 175, height = 120, units = "mm", dpi = 600, device = ragg::agg_tiff)

# Global weighted metrics
weighted_metrics <- function(obs, pred, w) {
  ok <- is.finite(obs) & is.finite(pred) & is.finite(w) & w > 0
  obs <- obs[ok]; pred <- pred[ok]; w <- w[ok]
  mu_w <- weighted.mean(obs, w)
  wRMSE <- sqrt(weighted.mean((obs - pred)^2, w))
  SS_res <- sum(w * (obs - pred)^2)
  SS_tot <- sum(w * (obs - mu_w)^2)
  wR2 <- 1 - SS_res / SS_tot
  tibble(wRMSE = wRMSE, wR2 = wR2)
}

glob_stats <- weighted_metrics(comp_plot_df$p_obs, comp_plot_df$p_hat, comp_plot_df$p_obs) %>%
  mutate(label = sprintf("wRMSE = %.3f\nwR² = %.2f", wRMSE, wR2), x = Inf, y = Inf)

p_global <- ggplot(comp_plot_df, aes(x = p_hat, y = p_obs, colour = Stress, shape = Diversity)) +
  geom_abline(linetype = 2) +
  geom_errorbarh(aes(xmin = p_lwr, xmax = p_upr), alpha = 0.3, height = 0, linewidth = 0.3) +
  geom_point(alpha = 0.8, size = 1.7) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Predicted relative abundance (posterior mean ± 90% CI)", y = "Observed relative abundance") +
  theme_classic(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 16, family = "Helvetica"), axis.text = element_text(size = 14, family = "Helvetica"), legend.title = element_text(size = 14, face = "bold", family = "Helvetica"), legend.text  = element_text(size = 12, family = "Helvetica"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)) +
  geom_text(data = glob_stats, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 12 / 2.845, family = "Helvetica")

ggsave("FIG_BAYES_comp_global_scatter_weighted.tiff", p_global, width = 170, height = 130, units = "mm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")

# Baseline bars
ggsave("FIG_BAYES_comp_baseline_bars.pdf", {
  comp_eq <- comp_plot_df %>% group_by(SampleID) %>% mutate(p_eq = 1 / n()) %>% ungroup()
  base_tbl <- tibble::tibble(
    model = c("Bayes softmax", "Equal-abundance"),
    RMSE  = c(sqrt(mean((comp_plot_df$p_obs - comp_plot_df$p_hat)^2, na.rm = TRUE)), sqrt(mean((comp_eq$p_obs - comp_eq$p_eq)^2, na.rm = TRUE)))
  )
  js_soft <- comp_plot_df %>% arrange(SampleID, Id) %>% group_by(SampleID) %>% summarise(JS = philentropy::JSD(rbind(p_obs, p_hat)), .groups = "drop")
  js_eq   <- comp_eq       %>% arrange(SampleID, Id) %>% group_by(SampleID) %>% summarise(JS = philentropy::JSD(rbind(p_obs, p_eq)),  .groups = "drop")
  base_tbl$JS <- c(mean(js_soft$JS, na.rm = TRUE), mean(js_eq$JS, na.rm = TRUE))
  base_long <- tidyr::pivot_longer(base_tbl, c(RMSE, JS), names_to = "metric", values_to = "value")
  ggplot(base_long, aes(metric, value, fill = model)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.55) +
    scale_fill_manual(values = c("Bayes softmax" = "grey20", "Equal-abundance" = "grey70")) +
    labs(x = NULL, y = "Error / divergence (lower = better)", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position = "top", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))
}, width = 110, height = 75, units = "mm")

# ==========================================================
# OD predictions: observed AWM vs Bayesian AWM vs baselines
# ==========================================================
message("[OD] Preparing community OD data ...")

comm_od <- read.csv("CommunityOD_All_Jan23_edit_zeros.csv") %>%
  mutate(
    Stress   = recode_stress(Stress),
    Com_Id   = paste0("Com", gsub("\\s+","", Community)),
    Diversity= as.character(Diversity_Level)
  ) %>%
  filter(Evo_Treatment != "Evo", t == 2) %>%
  transmute(
    Com_Id = gsub("\\s+","", Com_Id),
    Stress,
    Diversity,
    OD = as.numeric(OD)
  ) %>%
  mutate(Stress = factor(Stress), Diversity = factor(Diversity)) %>%
  droplevels()

lvl_stress <- levels(comm_od$Stress)
lvl_div    <- levels(comm_od$Diversity)

# Align ord to OD factor levels for safe joins
ord_od <- ord %>% mutate(Stress=factor(Stress, levels=lvl_stress), Diversity=factor(as.character(Diversity), levels=lvl_div))

# (1) Observed AWM (oracle ceiling)
message("[OD] Computing observed AWM (oracle ceiling) ...")
AWM_obs <- ord_od %>% group_by(SampleID, Com_Id, Stress, Diversity) %>% summarise(AWM_obs = sum(Abundance * g, na.rm = TRUE), .groups = "drop") %>% group_by(Com_Id, Stress, Diversity) %>% summarise(AWM_obs = mean(AWM_obs, na.rm = TRUE), .groups = "drop") %>% mutate(Stress=factor(Stress, levels=lvl_stress), Diversity=factor(Diversity, levels=lvl_div))

od_obs <- comm_od %>% left_join(AWM_obs, by = c("Com_Id","Stress","Diversity")) %>% filter(is.finite(OD), is.finite(AWM_obs))

m_obs <- lm(OD ~ AWM_obs * Stress + Diversity, data = od_obs)
od_obs$pred <- predict(m_obs, newdata = od_obs)
rmse_obs <- sqrt(mean((od_obs$OD - od_obs$pred)^2))
r2_obs   <- 1 - sum((od_obs$OD - od_obs$pred)^2) / sum((od_obs$OD - mean(od_obs$OD))^2)

p_obs <- ggplot(od_obs, aes(pred, OD, colour = Stress)) + geom_abline(linetype = 2) + geom_point(alpha = 0.9) + annotate("text", x = Inf, y = Inf, label = sprintf("RMSE = %.2f\nR² = %.2f", rmse_obs, r2_obs), hjust = 1.05, vjust = 1.4, size = 4) + labs(title="(1) Observed AWM × g", x = expression(Predicted~OD[600]), y = expression(OD[600])) + theme_bw()

# ---------- Helper: compute AWM per draw at Com×Stress×Diversity ----
compute_one_AWM_draw <- function(d){
  AWM_draw <- vector("list", N)
  for (n in seq_len(N)){
    seg <- idx_split[[n]]; s <- stress_id[n]; r <- if (USE_RICHNESS_KAPPA) rich_id[n] else 1L
    logits <- draws$kappa[d, s, r] * g_z[seg] + draws$delta[d, tax_of[seg]]
    p <- softmax(logits)
    AWM_draw[[n]] <- tibble(SampleID = SAMPLES$SampleID[n], AWM = sum(p * g_raw[seg]))
  }
  bind_rows(AWM_draw) %>% left_join(SAMPLES, by = "SampleID") %>% group_by(Com_Id, Stress, Diversity) %>% summarise(AWM = mean(AWM), .groups = "drop") %>% mutate(Stress=factor(Stress, levels=lvl_stress), Diversity=factor(as.character(Diversity), levels=lvl_div))
}

# Compute a point-estimate AWM (posterior mean across selected draws)
awm_point <- {
  tmp <- lapply(sel_mi, compute_one_AWM_draw)
  bind_rows(tmp, .id = "draw") %>% group_by(Com_Id, Stress, Diversity) %>% summarise(AWM_mean = mean(AWM), .groups = "drop")
}

# ---------------------- Nested-CV (optional) -------------------------
choose_od_formula <- function(df, awm_col = "AWM_mean", k_outer = 5, k_inner = 3) {
  set.seed(123)
  formulas <- list(
    f1 = as.formula(paste("OD ~", awm_col, "* Stress + Diversity")),
    f2 = as.formula(paste("OD ~", awm_col, "+ Stress + Diversity")),
    f3 = as.formula(paste("OD ~", awm_col)),
    f4 = as.formula("OD ~ Stress + Diversity")
  )
  # Create folds
  n <- nrow(df); idx <- sample.int(n)
  folds <- split(idx, cut(seq_along(idx), breaks = k_outer, labels = FALSE))
  mean_rmse <- rep(0, length(formulas))
  for (k in seq_along(folds)){
    test_idx  <- folds[[k]]
    train_idx <- setdiff(idx, test_idx)
    dtrain <- df[train_idx, , drop=FALSE]
    # inner CV to pick formula
    inner_idx <- sample(dtrain %>% rownames() %>% as.integer())
    inner_folds <- split(inner_idx, cut(seq_along(inner_idx), breaks = k_inner, labels = FALSE))
    rmse_mat <- matrix(NA_real_, nrow = length(formulas), ncol = length(inner_folds))
    for (j in seq_along(formulas)){
      for (ii in seq_along(inner_folds)){
        v_idx <- inner_folds[[ii]]; tr_idx <- setdiff(inner_idx, v_idx)
        fit <- tryCatch(lm(formulas[[j]], data = dtrain[tr_idx, , drop=FALSE]), error=function(e) NULL)
        if (is.null(fit)) next
        pred <- predict(fit, newdata = dtrain[v_idx, , drop=FALSE])
        rmse_mat[j, ii] <- sqrt(mean((dtrain$OD[v_idx] - pred)^2, na.rm=TRUE))
      }
    }
    inner_rmse <- rowMeans(rmse_mat, na.rm = TRUE)
    best_j <- which.min(inner_rmse)
    # train best on full train and eval on test
    fit_best <- lm(formulas[[best_j]], data = dtrain)
    pred_te  <- predict(fit_best, newdata = df[test_idx, , drop=FALSE])
    mean_rmse[best_j] <- mean_rmse[best_j] + sqrt(mean((df$OD[test_idx] - pred_te)^2, na.rm=TRUE))
  }
  best_overall <- which.max(mean_rmse == max(mean_rmse, na.rm = TRUE))
  formulas[[best_overall]]
}

# (2) Bayesian AWM via MI (trait-only)
message("[OD] Multiple-imputation over posterior compositions for AWM ...")

# Build base frame for formula selection using posterior-mean AWM
base_df <- comm_od %>% left_join(awm_point, by = c("Com_Id","Stress","Diversity")) %>% filter(is.finite(OD), is.finite(AWM_mean))

best_formula <- if (USE_NESTED_CV_FOR_OD) choose_od_formula(base_df, awm_col = "AWM_mean", k_outer = K_OUTER, k_inner = K_INNER) else as.formula("OD ~ AWM_mean * Stress + Diversity")

mi_preds <- vector("list", length(sel_mi)); keep <- 0L
for (ii in seq_along(sel_mi)){
  d <- sel_mi[ii]
  awm_d <- compute_one_AWM_draw(d)
  dat_d <- comm_od %>% left_join(awm_d, by = c("Com_Id","Stress","Diversity")) %>% filter(is.finite(OD), is.finite(AWM)) %>% rename(AWM_mean = AWM)
  if (nrow(dat_d) < 5) next
  m_d <- tryCatch(lm(best_formula, data = dat_d), error = function(e) NULL)
  if (is.null(m_d)) next
  keep <- keep + 1L
  mi_preds[[keep]] <- tibble(Com_Id=dat_d$Com_Id, Stress=dat_d$Stress, Diversity=dat_d$Diversity, OD=dat_d$OD, pred=as.numeric(predict(m_d, newdata = dat_d)))
}
if (keep == 0L) stop("[OD] MI failed: no usable fits. Check factor alignment and AWM availability.")
mi_preds <- mi_preds[seq_len(keep)] %>% bind_rows()
mi_summary <- mi_preds %>% group_by(Com_Id, Stress, Diversity, OD) %>% summarise(pred_mean = mean(pred), pred_sd = sd(pred), .groups = "drop")

rmse_bayes <- sqrt(mean((mi_summary$OD - mi_summary$pred_mean)^2))
r2_bayes   <- 1 - sum((mi_summary$OD - mi_summary$pred_mean)^2) / sum((mi_summary$OD - mean(mi_summary$OD))^2)

p_bayes <- ggplot(mi_summary, aes(pred_mean, OD, colour = Stress)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.9) +
  geom_errorbarh(aes(xmin = pred_mean - 1.96*pred_sd, xmax = pred_mean + 1.96*pred_sd), height = 0) +
  annotate("text", x = Inf, y = Inf, label = sprintf("RMSE = %.2f\nR² = %.2f", rmse_bayes, r2_bayes), hjust = 1.05, vjust = 1.4, size = 4) +
  labs(title="(2) Bayesian AWM × g (posterior MI)", x = expression(Predicted~OD[600]), y = expression(OD[600])) + theme_bw()

# --------------------- (3–5) Baselines (all × g) ---------------------
message("[OD] Building baselines (equal, random, stress-only) ...")

# Equal-abundance over planned members
eq_tbl <- design %>% left_join(all_growth, by = "Id", relationship = "many-to-many") %>% group_by(Com_Id, Stress) %>% summarise(AWM = mean(g, na.rm = TRUE), .groups="drop")

df2 <- comm_od %>% transmute(Com_Id, Stress, Diversity, OD) %>% left_join(eq_tbl, by=c("Com_Id","Stress")) %>% filter(is.finite(OD), is.finite(AWM))

best2 <- lm(OD ~ AWM + Stress + Diversity, data = df2)
df2$pred <- predict(best2, newdata = df2)
rmse2 <- sqrt(mean((df2$OD - df2$pred)^2))
r2_2  <- 1 - sum((df2$OD - df2$pred)^2)/sum((df2$OD - mean(df2$OD))^2)

p_eq <- ggplot(df2, aes(pred, OD, colour = Stress)) + geom_abline(linetype = 2) + geom_point(alpha = .9) + annotate("text", x = Inf, y = Inf, label = sprintf("RMSE = %.2f\nR² = %.2f", rmse2, r2_2), hjust = 1.05, vjust = 1.4, size = 4) + labs(title="(3) Equal-abundance baseline (×g)", x = expression(Predicted~OD[600]), y = expression(OD[600])) + theme_bw()

# Random composition: Dirichlet(1) weights
set.seed(42)
.rdir <- function(k){ x <- rgamma(k,1); x/sum(x) }

rand_tbl <- design %>% left_join(all_growth, by="Id", relationship = "many-to-many") %>% group_by(Com_Id, Stress) %>% summarise(AWM = { gvec <- g[is.finite(g)]; if (length(gvec)==0) NA_real_ else sum(.rdir(length(gvec)) * gvec) }, .groups="drop")

df3 <- comm_od %>% transmute(Com_Id, Stress, Diversity, OD) %>% left_join(rand_tbl, by=c("Com_Id","Stress")) %>% filter(is.finite(OD), is.finite(AWM))

best3 <- lm(OD ~ AWM + Stress + Diversity, data = df3)
df3$pred <- predict(best3, newdata = df3)
rmse3 <- sqrt(mean((df3$OD - df3$pred)^2))
r2_3  <- 1 - sum((df3$OD - df3$pred)^2)/sum((df3$OD - mean(df3$OD))^2)

p_rand <- ggplot(df3, aes(pred, OD, colour = Stress)) + geom_abline(linetype = 2) + geom_point(alpha = .9) + annotate("text", x = Inf, y = Inf, label = sprintf("RMSE = %.2f\nR² = %.2f", rmse3, r2_3), hjust = 1.05, vjust = 1.4, size = 4) + labs(title="(4) Random community baseline (×g)", x = expression(Predicted~OD[600]), y = expression(OD[600])) + theme_bw()

# Stress-only (no traits)
df4 <- comm_od %>% filter(is.finite(OD))
best4 <- lm(OD ~ Stress + Diversity, data = df4)
df4$pred <- predict(best4, newdata = df4)
rmse4 <- sqrt(mean((df4$OD - df4$pred)^2))
r2_4  <- 1 - sum((df4$OD - df4$pred)^2)/sum((df4$OD - mean(df4$OD))^2)

p_stress <- ggplot(df4, aes(pred, OD, colour = Stress)) + geom_abline(linetype = 2) + geom_point(alpha = .9) + annotate("text", x = Inf, y = Inf, label = sprintf("RMSE = %.2f\nR² = %.2f", rmse4, r2_4), hjust = 1.05, vjust = 1.4, size = 4) + labs(title="(5) Stress-only baseline", x = expression(Predicted~OD[600]), y = expression(OD[600])) + theme_bw()

# ===========================================================
# Apply Nature-style theme and save all AWM vs OD panels
# ===========================================================

theme_nat <- theme_classic(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 16, family = "Helvetica"), axis.text = element_text(size = 14, family = "Helvetica"), legend.title = element_text(size = 14, face = "bold", family = "Helvetica"), legend.text  = element_text(size = 12, family = "Helvetica"), plot.title = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))

plots <- list(FIG_OD_obs_AWM = p_obs, FIG_OD_bayes_AWM = p_bayes, FIG_OD_eq = p_eq, FIG_OD_rand = p_rand, FIG_OD_stress = p_stress)
for (nm in names(plots)) {
  p <- plots[[nm]] + theme_nat + theme(legend.position = "right")
  ggsave(paste0(nm, ".tiff"), p, width = 170, height = 130, units = "mm", dpi = 600, device = ragg::agg_tiff, compression = "lzw")
}

message("Done. ✔")



