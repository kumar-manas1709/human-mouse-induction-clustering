#!/usr/bin/env Rscript
#
# filter_orthologs.R
#
# Filters HumanRawCount.ods, HumanControlRawCount.ods, and MouseRawCount.ods
# down to only genes that are 1:1 human-mouse orthologs, using MGI's
# HOM_ProteinCoding.rpt ortholog report as the reference mapping.
#
# USAGE
# -----
# 1. Download the ortholog reference table from MGI:
#    https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt
#    Save it in the same folder as this script (or pass its path below).
#
# 2. Install dependencies (once):
#    install.packages(c("readODS", "dplyr", "readr", "stringr"))
#
# 3. Run:
#    Rscript filter_orthologs.R \
#        HumanRawCount.ods HumanControlRawCount.ods MouseRawCount.ods \
#        HOM_ProteinCoding.rpt
#
# OUTPUT
# ------
# - HumanRawCount_orthologs.csv
# - HumanControlRawCount_orthologs.csv
# - MouseRawCount_orthologs.csv            (mouse symbols kept as-is)
# - MouseRawCount_orthologs_humanized.csv  (mouse rows relabeled with their
#   human ortholog symbol, so you can align/merge/compare directly against
#   the human matrices row-for-row)
# - ortholog_mapping_used.csv              (the exact mouse<->human pairs applied)
#
# NOTES
# -----
# - HOM_ProteinCoding.rpt only reports strict 1:1 orthologs, 
#   genes with no listed ortholog, or with ambiguous (1:many / many:many) 
#   relationships, are dropped.
# - If you'd rather pull orthologs live from Ensembl instead of using the
#   static MGI file, see the biomaRt alternative at the bottom of this
#   script (commented out) -- useful if you want the very latest Ensembl
#   Compara calls rather than MGI's curated set.

suppressMessages({
  library(readODS)
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript filter_orthologs.R <human.ods> <human_control.ods> ",
    "<mouse.ods> <HOM_ProteinCoding.rpt>"
  )
}
human_path         <- args[1]
human_control_path <- args[2]
mouse_path         <- args[3]
ortholog_path      <- args[4]

# ---- helpers ---------------------------------------------------------

load_count_matrix <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- if (ext == "ods") {
    read_ods(path)
  } else if (ext %in% c("tsv", "txt")) {
    read_tsv(path, show_col_types = FALSE)
  } else {
    read_csv(path, show_col_types = FALSE)
  }
  names(df)[1] <- "gene_symbol"
  df$gene_symbol <- str_trim(as.character(df$gene_symbol))
  df
}

load_ortholog_table <- function(path) {
  raw <- read_tsv(path, show_col_types = FALSE)

  find_col <- function(...) {
    keywords <- tolower(c(...))
    hits <- names(raw)[
      sapply(tolower(names(raw)), function(n) all(sapply(keywords, function(k) str_detect(n, k))))
    ]
    if (length(hits) > 0) hits[1] else NA
  }

  mouse_col <- find_col("mouse", "symbol")
  human_col <- find_col("human", "symbol")

  if (is.na(mouse_col) || is.na(human_col)) {
    stop(
      "Could not find mouse/human symbol columns in the ortholog file.\n",
      "Columns found: ", paste(names(raw), collapse = ", "), "\n",
      "Rename the relevant columns to 'mouse_symbol' / 'human_symbol' and re-run."
    )
  }

  orth <- raw %>%
    select(mouse_symbol = all_of(mouse_col), human_symbol = all_of(human_col)) %>%
    filter(!is.na(mouse_symbol), !is.na(human_symbol)) %>%
    mutate(
      mouse_symbol = str_trim(mouse_symbol),
      human_symbol = str_trim(human_symbol)
    ) %>%
    distinct()

  # Keep strictly 1:1 pairs -- drop any symbol appearing more than once on
  # either side (ambiguous mappings).
  mouse_counts <- orth %>% count(mouse_symbol)
  human_counts <- orth %>% count(human_symbol)
  keep_mouse <- mouse_counts$mouse_symbol[mouse_counts$n == 1]
  keep_human <- human_counts$human_symbol[human_counts$n == 1]

  orth %>%
    filter(mouse_symbol %in% keep_mouse, human_symbol %in% keep_human)
}

# ---- main --------------------------------------------------------------

cat("Loading count matrices...\n")
human       <- load_count_matrix(human_path)
human_ctrl  <- load_count_matrix(human_control_path)
mouse       <- load_count_matrix(mouse_path)

cat("Loading ortholog mapping...\n")
orth <- load_ortholog_table(ortholog_path)
cat(sprintf("  %d clean 1:1 mouse<->human ortholog pairs\n", nrow(orth)))

human_orth_genes <- unique(orth$human_symbol)
mouse_orth_genes <- unique(orth$mouse_symbol)

human_f      <- human      %>% filter(gene_symbol %in% human_orth_genes)
human_ctrl_f <- human_ctrl %>% filter(gene_symbol %in% human_orth_genes)
mouse_f      <- mouse      %>% filter(gene_symbol %in% mouse_orth_genes)

cat(sprintf("Human:         %d -> %d genes with a mouse ortholog\n", nrow(human), nrow(human_f)))
cat(sprintf("Human control: %d -> %d genes with a mouse ortholog\n", nrow(human_ctrl), nrow(human_ctrl_f)))
cat(sprintf("Mouse:         %d -> %d genes with a human ortholog\n", nrow(mouse), nrow(mouse_f)))

write_csv(human_f,      "HumanRawCount_orthologs.csv")
write_csv(human_ctrl_f, "HumanControlRawCount_orthologs.csv")
write_csv(mouse_f,      "MouseRawCount_orthologs.csv")
write_csv(orth,         "ortholog_mapping_used.csv")

# "Humanized" mouse matrix: same data, rows relabeled with the human
# ortholog symbol so it lines up directly against the human matrices.
m2h <- setNames(orth$human_symbol, orth$mouse_symbol)
mouse_humanized <- mouse_f %>%
  mutate(gene_symbol = unname(m2h[gene_symbol]))
write_csv(mouse_humanized, "MouseRawCount_orthologs_humanized.csv")

cat("\nDone. Files written:\n")
cat("  HumanRawCount_orthologs.csv\n")
cat("  HumanControlRawCount_orthologs.csv\n")
cat("  MouseRawCount_orthologs.csv\n")
cat("  MouseRawCount_orthologs_humanized.csv\n")
cat("  ortholog_mapping_used.csv\n")

# -------------------------------------------------------------------------
# ALTERNATIVE: pull orthologs live from Ensembl via biomaRt instead of
# using the static MGI file. Uncomment and use in place of
# load_ortholog_table() above if you'd prefer this route (requires
# internet access and the biomaRt package: BiocManager::install("biomaRt")).
#
# library(biomaRt)
# human_mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
# mouse_mart <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl")
#
# orth_bm <- getLDS(
#   attributes  = c("hgnc_symbol"),
#   mart        = human_mart,
#   attributesL = c("mgi_symbol"),
#   martL       = mouse_mart,
#   uniqueRows  = TRUE
# )
# names(orth_bm) <- c("human_symbol", "mouse_symbol")
# orth_bm <- orth_bm %>% filter(human_symbol != "", mouse_symbol != "")
# -------------------------------------------------------------------------
