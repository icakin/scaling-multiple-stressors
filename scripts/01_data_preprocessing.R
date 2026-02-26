# ======================================================================
# Preprocessing ASV Data and Building Filtered Phyloseq Objects
#
# Overview:
# - Imports ASV counts, taxonomy, metadata, and a phylogenetic tree
# - Constructs a phyloseq object for the experimental dataset
# - Converts counts to relative abundance and filters low-abundance ASVs
# - Optionally merges two Pseudomonas ASVs into a single taxon
# - Subsets to Eight-species, non-evolved communities
# - Exports aligned ASV matrix, metadata, and taxonomy as RDS files
#
# Folder conventions:
# - Inputs (raw):      data/
# - Scripts:           scripts/
# - RDS outputs:       rds/        (via P_RDS())
#
# Outline:
# 1) Load helper functions and expected input files
# 2) Tidy taxonomy and rename ASV_10 genus (Pseudomonas2)
# 3) Import ASV counts and ensure correct orientation
# 4) Import and derive metadata (pH, Salinity; set factors)
# 5) Import phylogenetic tree and build phyloseq object
# 6) Transform to relative abundance and filter mean abundance > 0.5%
# 7) Optionally merge two Pseudomonas ASVs
# 8) Subset to Diversity == "Eight" and non-evo communities
# 9) Align ASV matrix and metadata by common samples
# 10) Save filtered phyloseq object, matrices, and taxonomy as RDS
# ======================================================================

# 01) Data preprocessing  -> saves aligned matrices & RDS
source("scripts/utils_functions.R"); ensure_packages()

# Expect in data/: ASV.tax.tsv, ASV.counts.tsv, Sequencing_metadata2.csv, ASV1-13_tree.txt

taxa_table <- read.csv(P_IN("ASV.tax.tsv"), sep = "\t")
taxa_table <- taxa_table %>%
  mutate(Genus = ifelse(row.names(.) == "ASV_10", "Pseudomonas2", Genus)) %>%
  microeco::tidy_taxonomy() %>% as.matrix()

OTU_data <- read.csv(P_IN("ASV.counts.tsv"), sep = "\t")
OTU_data <- as.matrix(t(OTU_data)) # samples as rows

sampleData <- read.csv(P_IN("Sequencing_metadata2.csv")) %>%
  mutate(
    pH       = ifelse(Stress %in% c("pH","pHSal","pHTemp","pHSalTemp"), 5.5, 7.2),
    Salinity = ifelse(Stress %in% c("Sal","pHSal","SalTemp","pHSalTemp"), 20, 0)
  )
sampleData$Temp     <- factor(sampleData$Temp)
sampleData$Salinity <- factor(sampleData$Salinity)
sampleData$pH       <- factor(sampleData$pH)
rownames(sampleData) <- sampleData$SampleID
sampleData <- dplyr::select(sampleData, -SampleID)

tree_data <- ape::read.tree(P_IN("ASV1-13_tree.txt"))

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

# Matrices & metadata (aligned)
ASV_t        <- as.data.frame(t(otu_table(ps_filter1)))
env_option_A <- as(phyloseq::sample_data(ps_filter1), "data.frame"); env_option_A$Id <- rownames(env_option_A)
common_samples     <- intersect(rownames(ASV_t), rownames(env_option_A))
ASV_t_clean        <- as.matrix(ASV_t[common_samples, , drop = FALSE]); storage.mode(ASV_t_clean) <- "double"
env_option_A_clean <- env_option_A[common_samples, , drop = FALSE]

# Save key objects
saveRDS(ps_filter1,       P_RDS("phyloseq_filtered.rds"))
saveRDS(ASV_t_clean,      P_RDS("ASV_matrix_clean.rds"))
saveRDS(env_option_A_clean, P_RDS("metadata_clean.rds"))
saveRDS(taxa_table,       P_RDS("taxa_table.rds"))

message("[01] Saved RDS: phyloseq_filtered.rds, ASV_matrix_clean.rds, metadata_clean.rds, taxa_table.rds")
