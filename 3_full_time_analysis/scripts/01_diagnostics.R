library(readr)
library(DESeq2)
library(Mfuzz)
library(Biobase)

source("pipeline_utils.R")  # provides out_path(dataset, filename)

# ============================================================================
# 01_diagnostics.R
#
# Part 1 of 2: everything up through cluster overlap evaluation --
# load -> clean -> DESeq2 normalize -> VST QC -> gene-level outlier check ->
# gene filtering -> overlap evaluation across candidate cluster counts.
#
# Does NOT run final clustering or produce cluster plots -- that's
# 02_clustering.R, which picks up from the "1_Output"/<dataset>/pipeline_state.rds
# files this script writes at the end. Splitting it this way means you can
# re-run clustering with a different c_final as many times as you like
# without ever re-running normalization/QC, which is the slow part.
#
# USAGE
# -----
#   Rscript 01_diagnostics.R
#
# (run from the folder containing pipeline_utils.R and your three matched
# CSVs -- output files are written under a folder named "1_Output",
# inside numbered per-dataset subfolders ("1_human", "2_mouse",
# "3_control" -- see pipeline_utils.R), created automatically if they
# don't exist)
#
# OUTPUT
# ------
# Per dataset, under "1_Output"/<numbered dataset folder>/ (see
# pipeline_utils.R for the exact folder names):
#   VST_Sample_Summary_Stats.csv
#   VST_Boxplot_By_Sample.pdf
#   VST_Mean_By_Sample.pdf
#   VST_Sample_Correlation_Matrix.csv
#   VST_Sample_Correlation_Heatmap.pdf
#   VST_Correlation_Boxplot_By_Sample.pdf
#   Gene_Mean_Expression_Histogram.pdf   (gene-level outlier check)
#   Gene_Mean_Expression_Outliers.csv    (genes >3 SD out on the log-scale mean)
#   Gene_CV_Histogram.pdf                (real CV distribution, for setting
#                                          cv_threshold from data, not a guess)
#   Mfuzz_Cluster_Overlap_Evaluation.pdf
#   pipeline_state.rds                   (for 02_clustering.R)
# ============================================================================

# ============================================================================
# CONFIG: the three matched files (same gene set, produced by
# match_gene_sets.R) that this script runs the pipeline on.
# ============================================================================

datasets <- list(
    human   = "HumanRawCount_matched.csv",
    mouse   = "MouseRawCount_matched.csv",
    control = "HumanControlRawCount_matched.csv"
)

# Samples flagged as outliers during VST QC (correlation boxplot / heatmap),
# to be dropped before normalization/filtering/clustering. Match substrings
# uniquely identifying the sample column -- no need for the exact full name.
outliers_to_exclude <- list(
    human   = c("NBplusNog_450", "NBplusNog_530", "NBplusNog_290", "NBplusNog_370"),
    mouse   = c("X.m_152", "X.m_468"),
    control = c("X.E8_410")
)

# ----------------------------------------------------------------------------
# Extracts the timepoint (minutes) from a column name.
#
# The timepoint is always the number immediately following the FIRST
# underscore in the name -- this was chosen deliberately, because some
# naming schemes have a digit embedded in the prefix itself before any
# real timepoint appears, e.g.:
#   "NBplusNog_360"            -> prefix "NBplusNog", first "_" -> 360
#   "X.m_580_rep_2."           -> prefix "X.m",        first "_" -> 580
#                                  (ignores the later "_2" replicate tag)
#   "X.E8_150mins."            -> prefix "X.E8" (contains digit "8", but
#                                  it's NOT preceded by "_"), first "_" -> 150
#   "X.E8_0._045220_0039683"   -> first "_" -> 0
#                                  (ignores the extra "_045220_0039683" junk
#                                  that got appended to this one column)
#
# If a column name has no "_<digits>" pattern at all, NA is returned and the
# calling code will warn you which column failed.
# ----------------------------------------------------------------------------
extract_time <- function(x) {
    m <- regmatches(x, regexpr("_[0-9]+", x))
    if (length(m) == 0 || m == "") return(NA_real_)
    as.numeric(sub("_", "", m))
}

# ============================================================================
# STAGE 1: Load the CSV and clean it into a numeric count matrix.
#
# `exclude_samples` is an optional character vector of substrings identifying
# columns to drop entirely (e.g. c("NBplusNog_450", "NBplusNog_530")) --
# used to remove samples flagged as outliers during VST QC, before they can
# influence normalization, filtering, or clustering.
# ============================================================================
load_and_clean <- function(csv_path, label, exclude_samples = character(0)) {
    count_data <- read_csv(csv_path, show_col_types = FALSE)
    count_data <- as.data.frame(count_data)
    rownames(count_data) <- count_data[[1]]
    count_data[[1]] <- NULL
    cat(label, "dimensions:", dim(count_data), "\n")

    if (length(exclude_samples) > 0) {
        drop_pattern <- paste(exclude_samples, collapse = "|")
        drop_cols <- grepl(drop_pattern, colnames(count_data))
        cat(
            label, "excluding", sum(drop_cols), "outlier sample(s):",
            paste(colnames(count_data)[drop_cols], collapse = ", "), "\n"
        )
        if (sum(drop_cols) != length(exclude_samples)) {
            warning(
                label, ": expected to exclude ", length(exclude_samples),
                " sample(s) but matched ", sum(drop_cols),
                " column(s) -- double check exclude_samples for typos."
            )
        }
        count_data <- count_data[, !drop_cols, drop = FALSE]
    }

    count_data[is.na(count_data)] <- 0
    valid_rows <- !is.na(rownames(count_data)) & rownames(count_data) != ""
    count_data <- count_data[valid_rows, ]

    count_matrix <- as.matrix(count_data)
    count_matrix <- round(count_matrix)
    count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
    cat(label, "genes with nonzero counts:", nrow(count_matrix), "\n")
    cat(label, "samples remaining:", ncol(count_matrix), "\n")

    count_matrix
}

# ============================================================================
# STAGE 2: Run DESeq2 size-factor normalization and VST.
# Returns a list with the DESeq2 object and the VST matrix.
# ============================================================================
normalize_deseq2 <- function(count_matrix) {
    col_data <- data.frame(
        condition = rep("sample", ncol(count_matrix)),
        row.names = colnames(count_matrix)
    )

    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = col_data,
        design = ~ 1
    )
    dds <- estimateSizeFactors(dds)

    vsd <- vst(dds, blind = TRUE)
    vst_matrix <- assay(vsd)

    list(dds = dds, vst_matrix = vst_matrix)
}

# ============================================================================
# STAGE 3: QC -- checks VST comparability across samples, ordered and
# colored by extracted timepoint. Writes CSV/PDF files under
# "1_Output"/<dataset>/. Returns the time-ordered VST matrix (handy downstream).
# ============================================================================
run_vst_qc <- function(vst_matrix, label) {
    cat("\n--- ", label, ": VST Sample Comparability QC ---\n")

    sample_time <- vapply(colnames(vst_matrix), extract_time, numeric(1))
    if (any(is.na(sample_time))) {
        warning(
            label, ": could not extract a time value from these column names: ",
            paste(colnames(vst_matrix)[is.na(sample_time)], collapse = ", "),
            " -- check the naming pattern. Falling back to original column order for those."
        )
        sample_time[is.na(sample_time)] <- seq_along(sample_time)[is.na(sample_time)]
    }

    time_order <- order(sample_time)
    vst_ordered <- vst_matrix[, time_order]
    time_ordered <- sample_time[time_order]

    time_colors <- colorRampPalette(c("blue", "red"))(length(unique(time_ordered)))[
        as.integer(factor(time_ordered))
    ]

    sample_stats <- data.frame(
        Dataset = label,
        sample = colnames(vst_ordered),
        time_minutes = time_ordered,
        mean   = colMeans(vst_ordered),
        median = apply(vst_ordered, 2, median),
        sd     = apply(vst_ordered, 2, sd)
    )
    print(sample_stats)
    write_csv(sample_stats, out_path(label, "VST_Sample_Summary_Stats.csv"))

    pdf(out_path(label, "VST_Boxplot_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.25), height = 6)
    boxplot(vst_ordered,
            las = 2,
            main = paste(label, ": VST-Normalized Expression by Sample (ordered by time)"),
            ylab = "VST value",
            col = time_colors,
            cex.axis = 0.6)
    abline(h = median(vst_ordered), col = "black", lty = 2)
    dev.off()

    pdf(out_path(label, "VST_Mean_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.25), height = 6)
    plot(time_ordered, sample_stats$mean,
         type = "b", pch = 19, col = time_colors,
         main = paste(label, ": Mean VST Value Over Time"),
         xlab = "Time (minutes)", ylab = "Mean VST value")
    abline(h = mean(sample_stats$mean), col = "grey40", lty = 2)
    dev.off()

    sample_cor <- cor(vst_ordered, method = "pearson")
    write.csv(sample_cor, file = out_path(label, "VST_Sample_Correlation_Matrix.csv"))

    pdf(out_path(label, "VST_Sample_Correlation_Heatmap.pdf"), width = 7, height = 7)
    heatmap(sample_cor,
            symm = TRUE,
            Rowv = NA, Colv = NA,
            main = paste(label, ": Sample-Sample Correlation (VST, ordered by time)"),
            margins = c(8, 8))
    dev.off()

    # Per-sample correlation boxplot -- the boxplot equivalent of the
    # heatmap above: for each sample, one box showing its correlation with
    # every OTHER sample (self-correlation, always 1, excluded). A box
    # sitting clearly below the rest is the same signal as a pale/white
    # row-or-column in the heatmap, just easier to read off precisely.
    # (Previously a separate script, correlation_boxplots.R -- folded in
    # here so it's never missed.)
    n_samples <- ncol(sample_cor)
    cor_list <- lapply(seq_len(n_samples), function(i) sample_cor[-i, i])
    names(cor_list) <- colnames(sample_cor)

    pdf(out_path(label, "VST_Correlation_Boxplot_By_Sample.pdf"),
        width = max(7, n_samples * 0.25), height = 6)
    boxplot(cor_list,
            las = 2,
            main = paste(label, ": Per-Sample Correlation to All Other Samples"),
            ylab = "Correlation with other samples",
            col = time_colors,
            cex.axis = 0.6)
    abline(h = median(unlist(cor_list)), col = "black", lty = 2)
    dev.off()

    cat(label, "QC files written to 1_Output/", dataset_folder_names[[label]], "/\n--- End QC ---\n\n", sep = "")

    vst_ordered
}

# ============================================================================
# STAGE 2b: Histogram of per-gene mean expression -- a GENE-level check
# (as opposed to Stage 3's SAMPLE-level checks). Each bar is one gene's
# mean normalized count across all samples. Since raw/normalized counts are
# heavily right-skewed (a handful of genes can be expressed 100-1000x more
# than the bulk), this is plotted on a log10 x-axis -- otherwise almost
# every gene piles into one bar at the low end and you can't see anything.
#
# What to look for:
#   - A single smooth-ish distribution (possibly bimodal: low-expressed
#     genes + a bump of highly-expressed ones) = normal/expected.
#   - A handful of genes isolated far out in the right tail, well separated
#     from the rest = worth a closer look (could be genuinely biologically
#     important, or could be a normalization artifact/contaminant).
#   - A spike of genes at exactly 0 = genes with no expression in this
#     dataset at all (these get removed by the mean>1 filter regardless).
#
# Writes a histogram PDF, and a CSV listing genes beyond `sd_threshold`
# standard deviations from the mean on the log10 scale -- a starting list
# of gene-level outlier candidates, analogous to the sample-level outliers
# flagged earlier from the correlation boxplots.
# ============================================================================
plot_mean_expression_histogram <- function(dds, label, sd_threshold = 3) {
    norm_counts <- counts(dds, normalized = TRUE)
    gene_means <- rowMeans(norm_counts)

    # log10(mean + 1) -- the "+1" avoids log(0) for unexpressed genes.
    log_means <- log10(gene_means + 1)

    pdf(out_path(label, "Gene_Mean_Expression_Histogram.pdf"), width = 8, height = 6)
    hist(log_means,
         breaks = 60,
         col = "steelblue",
         border = "white",
         main = paste(label, ": Per-Gene Mean Expression"),
         xlab = "log10(mean normalized count + 1)")
    abline(v = mean(log_means), col = "black", lty = 2)
    abline(v = mean(log_means) + sd_threshold * sd(log_means), col = "red", lty = 2)
    abline(v = mean(log_means) - sd_threshold * sd(log_means), col = "red", lty = 2)
    legend("topright", legend = c("mean", paste0("\u00b1", sd_threshold, " SD")),
           col = c("black", "red"), lty = 2, bty = "n")
    dev.off()
    cat(label, "wrote Gene_Mean_Expression_Histogram.pdf\n")

    # Flag genes far out in either tail on the log scale.
    z <- (log_means - mean(log_means)) / sd(log_means)
    outlier_genes <- data.frame(
        Dataset = label,
        Gene = names(gene_means),
        mean_normalized_count = gene_means,
        log10_mean = log_means,
        z_score = z
    )
    outlier_genes <- outlier_genes[abs(outlier_genes$z_score) > sd_threshold, ]
    outlier_genes <- outlier_genes[order(-abs(outlier_genes$z_score)), ]

    write_csv(outlier_genes, out_path(label, "Gene_Mean_Expression_Outliers.csv"))
    cat(
        label, ":", nrow(outlier_genes),
        "gene(s) beyond", sd_threshold, "SD on log-scale mean expression -- see",
        "Gene_Mean_Expression_Outliers.csv\n"
    )

    invisible(outlier_genes)
}

# ============================================================================
# STAGE 4: Filter genes by mean expression, then by coefficient of
# variation -- CV threshold is now the 75th percentile of each dataset's
# OWN CV distribution (see Gene_CV_Histogram.pdf), not a fixed guessed
# number, since CV scale differs meaningfully between datasets. This
# guarantees roughly the top 25% most variable mean>1 genes are kept, so
# the top-2000-by-variance fallback below should no longer trigger in
# practice -- it's kept only as a safety net.
# ============================================================================
filter_genes <- function(dds, vst_matrix, label) {
    norm_counts <- counts(dds, normalized = TRUE)
    keep_mean <- rowMeans(norm_counts) > 1
    cat(label, "genes with mean > 1:", sum(keep_mean), "\n")

    # Coefficient of variation, computed on VST values. NOTE: VST is a
    # log-like transform that deliberately compresses the scale and flattens
    # the mean-variance relationship -- so CV values here are naturally much
    # smaller than you'd see on raw/linear counts. A fixed guess like 0.4
    # can end up unreachable regardless of data quality. To avoid guessing,
    # print the actual CV distribution (quantiles + a histogram) for genes
    # that already passed the mean filter, so the threshold below can be
    # set from what the data actually looks like.
    cv_all <- apply(vst_matrix, 1, sd) / abs(rowMeans(vst_matrix))
    cv_mean_passed <- cv_all[keep_mean]

    cat(label, "CV distribution (genes with mean > 1):\n")
    print(quantile(cv_mean_passed, probs = c(0, .25, .5, .75, .9, .95, .99, 1)))

    pdf(out_path(label, "Gene_CV_Histogram.pdf"), width = 8, height = 6)
    hist(cv_mean_passed,
         breaks = 60,
         col = "darkolivegreen",
         border = "white",
         main = paste(label, ": Per-Gene CV (mean > 1 genes only)"),
         xlab = "Coefficient of variation (sd / mean, on VST scale)")
    abline(v = quantile(cv_mean_passed, 0.75), col = "red", lty = 2)
    legend("topright", legend = "75th percentile", col = "red", lty = 2, bty = "n")
    dev.off()
    cat(label, "wrote Gene_CV_Histogram.pdf\n")

    # Percentile-based threshold instead of a fixed absolute number: CV
    # values on the VST scale differ meaningfully between datasets (compare
    # the printed quantiles above -- mouse tops out around 0.15, human and
    # control around 0.07), so one fixed cutoff can't be right for all
    # three at once. cv_percentile = 0.75 keeps the most variable quarter
    # of mean>1 genes for THIS dataset, whatever CV value that happens to
    # correspond to -- calibrated from the data instead of guessed.
    cv_percentile <- 0.75
    cv_threshold <- quantile(cv_mean_passed, cv_percentile, names = FALSE)
    keep_cv <- cv_all > cv_threshold
    cat(label, "genes with CV >", round(cv_threshold, 4),
        paste0("(", cv_percentile * 100, "th percentile)"), ":", sum(keep_cv), "\n")

    final_vst_data <- vst_matrix[keep_mean & keep_cv, ]
    cat(label, "genes passing BOTH filters:", nrow(final_vst_data), "\n")

    if (nrow(final_vst_data) == 0) {
        warning(label, ": CV filter removed all genes! Falling back to top 2,000 variable genes (with mean > 1).")
        vst_high_mean <- vst_matrix[keep_mean, ]
        gene_var <- apply(vst_high_mean, 1, var)
        top_n <- min(2000, nrow(vst_high_mean))
        keep_top <- order(gene_var, decreasing = TRUE)[1:top_n]
        final_vst_data <- vst_high_mean[keep_top, ]
        cat(label, "fallback: retained", nrow(final_vst_data), "genes by variance.\n")
    }

    final_vst_data
}

# ============================================================================
# STAGE 5a: Evaluate cluster overlap across candidate c values (2/3/4/6/8/10),
# writing one diagnostic PDF (sensitive threshold) to help you pick c_final
# by eye. Returns m_opt (the fuzzifier estimated by Mfuzz for this dataset).
# ============================================================================
evaluate_cluster_overlap <- function(eset_r, label) {
    m_opt <- mestimate(eset_r)
    cat(label, "m_opt:", m_opt, "\n")

    # Only the sensitive (thres = 0.1) overlap evaluation is generated --
    # thres = 0.5 is a much higher correlation bar to clear, so it tends to
    # under-flag real overlap between clusters.
    pdf(out_path(label, "Mfuzz_Cluster_Overlap_Evaluation.pdf"), width = 7, height = 7)
    for (k in c(2, 3, 4, 6, 8, 10)) {
        set.seed(42)
        cl_test <- mfuzz(eset_r, c = k, m = m_opt)
        O <- overlap(cl_test)
        overlap.plot(cl_test, over = O, thres = 0.1)
        title(main = paste(label, ": Cluster Overlap for c =", k))
    }
    dev.off()

    m_opt
}

# ============================================================================
# ORCHESTRATOR: runs stages 1-5a, in order, for one dataset, then saves
# everything 02_clustering.R will need to a single RDS file.
# ============================================================================
run_diagnostics <- function(csv_path, label, exclude_samples = character(0)) {
    cat("\n================ ", label, " ================\n")

    count_matrix <- load_and_clean(csv_path, label, exclude_samples = exclude_samples)

    norm <- normalize_deseq2(count_matrix)
    dds <- norm$dds
    vst_matrix <- norm$vst_matrix

    run_vst_qc(vst_matrix, label)

    plot_mean_expression_histogram(dds, label)

    final_vst_data <- filter_genes(dds, vst_matrix, label)

    eset <- ExpressionSet(assayData = as.matrix(final_vst_data))
    eset_r <- standardise(eset)

    m_opt <- evaluate_cluster_overlap(eset_r, label)

    state_path <- out_path(label, "pipeline_state.rds")
    saveRDS(
        list(eset_r = eset_r, final_vst_data = final_vst_data, m_opt = m_opt, label = label),
        state_path
    )
    cat(label, ": diagnostics complete. Saved state to", state_path, "\n")

    invisible(NULL)
}

# ============================================================================
# Run diagnostics once per dataset
# ============================================================================

for (label in names(datasets)) {
    run_diagnostics(
        datasets[[label]], label,
        exclude_samples = outliers_to_exclude[[label]]
    )
}

cat("\nAll datasets diagnosed:", paste(names(datasets), collapse = ", "), "\n")
cat("Files are under \"1_Output\"/<dataset>/. Next: review each dataset's\n")
cat("Mfuzz_Cluster_Overlap_Evaluation.pdf, choose a c_final, then run 02_clustering.R.\n")
