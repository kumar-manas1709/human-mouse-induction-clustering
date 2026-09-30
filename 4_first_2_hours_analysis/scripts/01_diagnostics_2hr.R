library(readr)
library(DESeq2)
library(Mfuzz)
library(Biobase)

# ============================================================================
# 01_diagnostics_2hr.R
#
# Variant of ../../scripts/01_diagnostics.R, restricted to samples at or
# before MAX_TIME_MINUTES (120 = first two hours) for all three datasets.
# Deliberately self-contained (does not source the main pipeline_utils.R)
# so this exploratory analysis can never accidentally affect the main
# pipeline's output folder or vice versa.
#
# Same stages, same logic, same outlier exclusions as the main pipeline --
# the only functional difference is the time-window filter applied right
# after loading, and a separate output folder so results never mix with
# the full-timecourse run.
#
# CAVEAT: at <=120 minutes, human and control have only 13 timepoints each
# (10-minute intervals) and mouse has 31 (4-minute intervals) -- a much
# smaller sample than the full run. Treat this as exploratory; statistical
# power (e.g. the CV-percentile gene filter, cluster stability) is reduced
# accordingly.
#
# USAGE
# -----
#   Rscript 01_diagnostics_2hr.R
#
# (run from this folder -- first_2_hours_analysis/scripts/ -- so the
# relative path to ../../data/matched/ resolves correctly)
#
# OUTPUT
# ------
# Per dataset, under 1_Output_First2Hours/<dataset>/ (relative to wherever
# this script is run from):
#   VST_Sample_Summary_Stats.csv, VST_Boxplot_By_Sample.pdf,
#   VST_Mean_By_Sample.pdf, VST_Sample_Correlation_Matrix.csv,
#   VST_Sample_Correlation_Heatmap.pdf, VST_Correlation_Boxplot_By_Sample.pdf,
#   Gene_Mean_Expression_Histogram.pdf, Gene_Mean_Expression_Outliers.csv,
#   Gene_CV_Histogram.pdf, Mfuzz_Cluster_Overlap_Evaluation.pdf,
#   pipeline_state.rds (for 02_clustering_2hr.R)
# ============================================================================

MAX_TIME_MINUTES <- 120

# ============================================================================
# CONFIG: same matched files as the main pipeline, referenced relative to
# this script's location.
# ============================================================================

datasets <- list(
    human   = "../../data/matched/HumanRawCount_matched.csv",
    mouse   = "../../data/matched/MouseRawCount_matched.csv",
    control = "../../data/matched/HumanControlRawCount_matched.csv"
)

# Same exclusions as the main pipeline (applied for consistency -- none of
# these happen to fall within the <=120min window for this dataset, but the
# logic is kept identical so this subset is filtered by the same standard).
outliers_to_exclude <- list(
    human   = c("NBplusNog_450", "NBplusNog_530", "NBplusNog_290", "NBplusNog_370"),
    mouse   = c("X.m_152", "X.m_468"),
    control = c("X.E8_410")
)

dataset_folder_names <- c(human = "1_human", mouse = "2_mouse", control = "3_control")

out_path <- function(dataset, filename) {
    folder_name <- if (dataset %in% names(dataset_folder_names)) dataset_folder_names[[dataset]] else dataset
    dir <- file.path("1_Output_First2Hours", folder_name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file.path(dir, filename)
}

# ----------------------------------------------------------------------------
# Same timepoint-extraction rule as the main pipeline: the number right
# after the FIRST underscore in the column name.
# ----------------------------------------------------------------------------
extract_time <- function(x) {
    m <- regmatches(x, regexpr("_[0-9]+", x))
    if (length(m) == 0 || m == "") return(NA_real_)
    as.numeric(sub("_", "", m))
}

# ============================================================================
# STAGE 1: Load, clean, exclude known outlier samples, THEN restrict to
# samples at or before MAX_TIME_MINUTES.
# ============================================================================
load_and_clean <- function(csv_path, label, exclude_samples = character(0),
                            max_time = MAX_TIME_MINUTES) {
    count_data <- read_csv(csv_path, show_col_types = FALSE)
    count_data <- as.data.frame(count_data)
    rownames(count_data) <- count_data[[1]]
    count_data[[1]] <- NULL
    cat(label, "dimensions (before time filter):", dim(count_data), "\n")

    if (length(exclude_samples) > 0) {
        drop_pattern <- paste(exclude_samples, collapse = "|")
        drop_cols <- grepl(drop_pattern, colnames(count_data))
        cat(label, "excluding", sum(drop_cols), "outlier sample(s):",
            paste(colnames(count_data)[drop_cols], collapse = ", "), "\n")
        count_data <- count_data[, !drop_cols, drop = FALSE]
    }

    sample_time <- vapply(colnames(count_data), extract_time, numeric(1))
    if (any(is.na(sample_time))) {
        warning(label, ": could not extract a time value from: ",
                paste(colnames(count_data)[is.na(sample_time)], collapse = ", "))
    }
    keep_time <- !is.na(sample_time) & sample_time <= max_time
    cat(label, "samples at <=", max_time, "minutes:", sum(keep_time),
        "(dropped", sum(!keep_time), "later timepoints)\n")
    count_data <- count_data[, keep_time, drop = FALSE]

    count_data[is.na(count_data)] <- 0
    valid_rows <- !is.na(rownames(count_data)) & rownames(count_data) != ""
    count_data <- count_data[valid_rows, ]

    count_matrix <- as.matrix(count_data)
    count_matrix <- round(count_matrix)
    count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
    cat(label, "genes with nonzero counts:", nrow(count_matrix), "\n")
    cat(label, "final sample count:", ncol(count_matrix), "\n")

    count_matrix
}

# ============================================================================
# STAGE 2: DESeq2 normalization + VST (identical to the main pipeline;
# re-run on the subset so size factors reflect only these samples).
# ============================================================================
normalize_deseq2 <- function(count_matrix) {
    col_data <- data.frame(condition = rep("sample", ncol(count_matrix)),
                            row.names = colnames(count_matrix))
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ 1)
    dds <- estimateSizeFactors(dds)
    vsd <- vst(dds, blind = TRUE)
    list(dds = dds, vst_matrix = assay(vsd))
}

# ============================================================================
# STAGE 3: Same VST QC as the main pipeline (boxplot, mean-over-time,
# correlation heatmap, correlation boxplot), on the subset.
# ============================================================================
run_vst_qc <- function(vst_matrix, label) {
    cat("\n--- ", label, ": VST Sample Comparability QC (<=", MAX_TIME_MINUTES, "min) ---\n")

    sample_time <- vapply(colnames(vst_matrix), extract_time, numeric(1))
    time_order <- order(sample_time)
    vst_ordered <- vst_matrix[, time_order]
    time_ordered <- sample_time[time_order]

    time_colors <- colorRampPalette(c("blue", "red"))(length(unique(time_ordered)))[
        as.integer(factor(time_ordered))
    ]

    sample_stats <- data.frame(
        Dataset = label, sample = colnames(vst_ordered), time_minutes = time_ordered,
        mean = colMeans(vst_ordered), median = apply(vst_ordered, 2, median),
        sd = apply(vst_ordered, 2, sd)
    )
    print(sample_stats)
    write_csv(sample_stats, out_path(label, "VST_Sample_Summary_Stats.csv"))

    pdf(out_path(label, "VST_Boxplot_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.3), height = 6)
    boxplot(vst_ordered, las = 2, main = paste(label, ": VST-Normalized Expression, <=", MAX_TIME_MINUTES, "min"),
            ylab = "VST value", col = time_colors, cex.axis = 0.7)
    abline(h = median(vst_ordered), col = "black", lty = 2)
    dev.off()

    pdf(out_path(label, "VST_Mean_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.3), height = 6)
    plot(time_ordered, sample_stats$mean, type = "b", pch = 19, col = time_colors,
         main = paste(label, ": Mean VST Value Over Time, <=", MAX_TIME_MINUTES, "min"),
         xlab = "Time (minutes)", ylab = "Mean VST value")
    abline(h = mean(sample_stats$mean), col = "grey40", lty = 2)
    dev.off()

    sample_cor <- cor(vst_ordered, method = "pearson")
    write.csv(sample_cor, file = out_path(label, "VST_Sample_Correlation_Matrix.csv"))

    pdf(out_path(label, "VST_Sample_Correlation_Heatmap.pdf"), width = 7, height = 7)
    heatmap(sample_cor, symm = TRUE, Rowv = NA, Colv = NA,
            main = paste(label, ": Sample-Sample Correlation, <=", MAX_TIME_MINUTES, "min"),
            margins = c(8, 8))
    dev.off()

    n_samples <- ncol(sample_cor)
    cor_list <- lapply(seq_len(n_samples), function(i) sample_cor[-i, i])
    names(cor_list) <- colnames(sample_cor)
    pdf(out_path(label, "VST_Correlation_Boxplot_By_Sample.pdf"), width = max(7, n_samples * 0.3), height = 6)
    boxplot(cor_list, las = 2, main = paste(label, ": Per-Sample Correlation, <=", MAX_TIME_MINUTES, "min"),
            ylab = "Correlation with other samples", col = time_colors, cex.axis = 0.7)
    abline(h = median(unlist(cor_list)), col = "black", lty = 2)
    dev.off()

    cat(label, "QC files written.\n--- End QC ---\n\n")
    vst_ordered
}

# ============================================================================
# STAGE 2b: Gene-level mean expression histogram (same logic as main).
# ============================================================================
plot_mean_expression_histogram <- function(dds, label, sd_threshold = 3) {
    norm_counts <- counts(dds, normalized = TRUE)
    gene_means <- rowMeans(norm_counts)
    log_means <- log10(gene_means + 1)

    pdf(out_path(label, "Gene_Mean_Expression_Histogram.pdf"), width = 8, height = 6)
    hist(log_means, breaks = 60, col = "steelblue", border = "white",
         main = paste(label, ": Per-Gene Mean Expression, <=", MAX_TIME_MINUTES, "min"),
         xlab = "log10(mean normalized count + 1)")
    abline(v = mean(log_means), col = "black", lty = 2)
    abline(v = mean(log_means) + sd_threshold * sd(log_means), col = "red", lty = 2)
    abline(v = mean(log_means) - sd_threshold * sd(log_means), col = "red", lty = 2)
    legend("topright", legend = c("mean", paste0("\u00b1", sd_threshold, " SD")), col = c("black", "red"), lty = 2, bty = "n")
    dev.off()

    z <- (log_means - mean(log_means)) / sd(log_means)
    outlier_genes <- data.frame(Dataset = label, Gene = names(gene_means),
                                 mean_normalized_count = gene_means, log10_mean = log_means, z_score = z)
    outlier_genes <- outlier_genes[abs(outlier_genes$z_score) > sd_threshold, ]
    outlier_genes <- outlier_genes[order(-abs(outlier_genes$z_score)), ]
    write_csv(outlier_genes, out_path(label, "Gene_Mean_Expression_Outliers.csv"))
    cat(label, ":", nrow(outlier_genes), "gene(s) beyond", sd_threshold, "SD\n")
    invisible(outlier_genes)
}

# ============================================================================
# STAGE 4: Same two-step gene filter as the main pipeline (mean>1, then
# top-quartile CV, threshold set from THIS subset's own distribution).
# ============================================================================
filter_genes <- function(dds, vst_matrix, label) {
    norm_counts <- counts(dds, normalized = TRUE)
    keep_mean <- rowMeans(norm_counts) > 1
    cat(label, "genes with mean > 1:", sum(keep_mean), "\n")

    cv_all <- apply(vst_matrix, 1, sd) / abs(rowMeans(vst_matrix))
    cv_mean_passed <- cv_all[keep_mean]

    cat(label, "CV distribution (genes with mean > 1):\n")
    print(quantile(cv_mean_passed, probs = c(0, .25, .5, .75, .9, .95, .99, 1)))

    pdf(out_path(label, "Gene_CV_Histogram.pdf"), width = 8, height = 6)
    hist(cv_mean_passed, breaks = 60, col = "darkolivegreen", border = "white",
         main = paste(label, ": Per-Gene CV, <=", MAX_TIME_MINUTES, "min (mean > 1 genes only)"),
         xlab = "Coefficient of variation (sd / mean, on VST scale)")
    abline(v = quantile(cv_mean_passed, 0.75), col = "red", lty = 2)
    legend("topright", legend = "75th percentile", col = "red", lty = 2, bty = "n")
    dev.off()

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
# STAGE 5a: Cluster overlap evaluation (thres = 0.1 only, same as main).
# ============================================================================
evaluate_cluster_overlap <- function(eset_r, label) {
    m_opt <- mestimate(eset_r)
    cat(label, "m_opt:", m_opt, "\n")

    pdf(out_path(label, "Mfuzz_Cluster_Overlap_Evaluation.pdf"), width = 7, height = 7)
    for (k in c(2, 3, 4, 6, 8, 10)) {
        set.seed(42)
        cl_test <- mfuzz(eset_r, c = k, m = m_opt)
        O <- overlap(cl_test)
        overlap.plot(cl_test, over = O, thres = 0.1)
        title(main = paste(label, ": Cluster Overlap for c =", k, "(<=", MAX_TIME_MINUTES, "min)"))
    }
    dev.off()
    m_opt
}

# ============================================================================
# ORCHESTRATOR
# ============================================================================
run_diagnostics <- function(csv_path, label, exclude_samples = character(0)) {
    cat("\n================ ", label, " (<=", MAX_TIME_MINUTES, "min) ================\n")
    count_matrix <- load_and_clean(csv_path, label, exclude_samples = exclude_samples)
    norm <- normalize_deseq2(count_matrix)
    run_vst_qc(norm$vst_matrix, label)
    plot_mean_expression_histogram(norm$dds, label)
    final_vst_data <- filter_genes(norm$dds, norm$vst_matrix, label)

    eset <- ExpressionSet(assayData = as.matrix(final_vst_data))
    eset_r <- standardise(eset)
    m_opt <- evaluate_cluster_overlap(eset_r, label)

    state_path <- out_path(label, "pipeline_state.rds")
    saveRDS(list(eset_r = eset_r, final_vst_data = final_vst_data, m_opt = m_opt, label = label), state_path)
    cat(label, ": diagnostics complete. Saved state to", state_path, "\n")
    invisible(NULL)
}

for (label in names(datasets)) {
    run_diagnostics(datasets[[label]], label, exclude_samples = outliers_to_exclude[[label]])
}

cat("\nAll datasets diagnosed (<=", MAX_TIME_MINUTES, "min). Files under 1_Output_First2Hours/<dataset>/.\n")
cat("Next: review overlap PDFs, choose c_final, then run 02_clustering_2hr.R.\n")
