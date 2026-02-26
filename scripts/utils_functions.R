# utils_functions.R
# Shared helpers: packages, paths, plotting, palettes, seeds, dir creation

# ---- Packages ----
ensure_packages <- function() {
  pkgs <- c(
    "magrittr","dplyr","tidyr","tibble","ggplot2","ggpubr","scales","ggrepel","patchwork",
    "ape","phyloseq","vegan","pairwiseAdonis","ggtree","stringr","picante","phangorn","microeco","reshape2",
    "broom","car","effectsize","emmeans","readr","permute","indicspecies","ragg",
    "rstan","loo","philentropy"
  )
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install)
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!"pairwiseAdonis" %in% installed.packages()[,1]) {
    devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

# ---- Paths ----
DIR_INPUT  <- "data"
DIR_TABLES <- "results/tables"
DIR_FIGS   <- "results/figures"
DIR_RDS    <- "results/rds"

invisible(lapply(
  c(DIR_INPUT, DIR_TABLES, DIR_FIGS, DIR_RDS),
  function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
))

P_IN  <- function(...) file.path(DIR_INPUT,  ...)
P_TAB <- function(...) file.path(DIR_TABLES, ...)
P_FIG <- function(...) file.path(DIR_FIGS,   ...)
P_RDS <- function(...) file.path(DIR_RDS,    ...)

# ---- Global settings ----
set.seed(123)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(), stringsAsFactors = FALSE)

# ---- Plot tag ----
label_plot <- function(p, tag) {
  p + labs(tag = tag) +
    theme(
      plot.tag = element_text(size = 18, face = "bold", family = "Helvetica"),
      plot.tag.position = c(0.01, 0.98)
    )
}

# ---- Stress labels & palettes ----
levels_short <<- c("Control","pH","Sal","Temp","pHSal","pHTemp","SalTemp","pHSalTemp")

pretty_map   <<- c(
  "Control"   = "Control",
  "pH"        = "pH",
  "Sal"       = "Salinity",
  "Temp"      = "Temperature",
  "pHSal"     = "pH x Salinity",
  "pHTemp"    = "pH x Temperature",
  "SalTemp"   = "Salinity x Temperature",
  "pHSalTemp" = "pH x Salinity x Temperature"
)

pretty_levels <<- unname(pretty_map[levels_short])

pal_short <<- c(
  "Control"="darkgrey",
  "pH"     ="#56B4E9",
  "Sal"    ="#F0E442",
  "Temp"   ="#FF7043",
  "pHSal"  ="#009E73",
  "pHTemp" ="#CC79A7",
  "SalTemp"="chocolate4",
  "pHSalTemp"="#3949AB"
)

pal_pretty <<- setNames(unname(pal_short), pretty_levels)

# ---- Misc helpers ----
trimf    <- function(x) trimws(as.character(x))
safe_num <- function(x) suppressWarnings(as.numeric(x))

# reorder_within helpers (used in Fig 2)
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  ggplot2::scale_x_discrete(
    labels = function(x) gsub(paste0(sep, ".*$"), "", x),
    ...
  )
}

# ---- Softmax helper (used in Fig 3 script) ----
softmax <- function(x) {
  x <- as.numeric(x)
  x <- x - max(x, na.rm = TRUE)  # numerical stability
  ex <- exp(x)
  ex / sum(ex, na.rm = TRUE)
}
