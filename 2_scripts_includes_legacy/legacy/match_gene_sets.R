
suppressMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript match_gene_sets.R <human_orthologs.csv> <mouse_orthologs_humanized.csv> [human_control_orthologs.csv]")
}
human_path   <- args[1]
mouse_path   <- args[2]
control_path <- if (length(args) >= 3) args[3] else NA

cat("Loading files...\n")
human <- read_csv(human_path, show_col_types = FALSE)
mouse <- read_csv(mouse_path, show_col_types = FALSE)
control <- if (!is.na(control_path)) read_csv(control_path, show_col_types = FALSE) else NULL

human <- human %>% distinct(gene_symbol, .keep_all = TRUE)
mouse <- mouse %>% distinct(gene_symbol, .keep_all = TRUE)
if (!is.null(control)) {
  control <- control %>% distinct(gene_symbol, .keep_all = TRUE)
}

common_genes <- intersect(human$gene_symbol, mouse$gene_symbol)
cat(sprintf(
  "Human: %d genes, Mouse: %d genes -> %d genes in common\n",
  nrow(human), nrow(mouse), length(common_genes)
))

human_matched <- human %>%
  filter(gene_symbol %in% common_genes) %>%
  arrange(gene_symbol)

mouse_matched <- mouse %>%
  filter(gene_symbol %in% common_genes) %>%
  arrange(gene_symbol)

# Sanity check: confirm the two files now line up gene-for-gene, row-for-row.
stopifnot(identical(human_matched$gene_symbol, mouse_matched$gene_symbol))

write_csv(human_matched, "HumanRawCount_matched.csv")
write_csv(mouse_matched, "MouseRawCount_matched.csv")


if (!is.null(control)) {
  control_matched <- control %>%
    filter(gene_symbol %in% common_genes) %>%
    arrange(gene_symbol)

  write_csv(control_matched, "HumanControlRawCount_matched.csv")


  if (!identical(control_matched$gene_symbol, human_matched$gene_symbol)) {
    cat(
      "\nNote: the control file's matched genes differ slightly from the ",
      "main human/mouse set -- this can happen if the control file is ",
      "missing a few of those genes entirely. Check for NAs/gaps before ",
      "using it in a row-for-row comparison.\n",
      sep = ""
    )
  }
}
