library(readr)
library(DESeq2)
library(Mfuzz)
library(Biobase)

# Routes every output file from this run to "1_Output_120min" instead of
# the main pipeline's "1_Output" -- must be set BEFORE sourcing
# pipeline_utils.R. This keeps this analysis fully separate from the
# full-timecourse pipeline; nothing here can overwrite or interfere with
# 01_diagnostics.R's output, even if both are run from the same folder.
OUTPUT_ROOT <- "1_Output_120min"
source("pipeline_utils.R")  # provides out_path(dataset, filename)

# ============================================================================
# 01_diagnostics_120min.R
#
# EARLY-WINDOW VARIANT of 01_diagnostics.R: identical pipeline logic
# (load -> clean -> DESeq2 normalize -> VST QC -> gene-level outlier check ->
# gene filtering -> overlap evaluation), but restricted to samples from the
# first 2 hours (<=120 minutes) of each time course only. This is a FULL
# separate run: DESeq2 normalization (size factors, VST) is recomputed
# from scratch on just this subset, not sliced from the full-timecourse
# results, since size factors depend on which samples are included.
#
# Applies to all three datasets (human, mouse, control).
#
# WHY THIS IS A SEPARATE SCRIPT rather than a flag on the original:
# keeping the two fully independent means the validated full-timecourse
# pipeline (01_diagnostics.R / 02_clustering.R) is never at risk of being
# broken by this addition. This script only reads the same input CSVs
# and pipeline_utils.R; it does not modify either.
#
# USAGE
# -----
#   Rscript 01_diagnostics_120min.R
#
# (run from the same folder as 01_diagnostics.R, pipeline_utils.R, and your
# three matched CSVs -- output files are written under "1_Output_120min",
# with the same per-dataset subfolder structure as the main pipeline)
#
# OUTPUT
# ------
# Same file set as 01_diagnostics.R, per dataset, under
# "1_Output_120min"/<numbered dataset folder>/ -- see that script's header
# for the full list. Console output additionally reports how many samples
# were dropped by the 120-minute cutoff, separately from outlier exclusion.
# ============================================================================

# ============================================================================
# CONFIG
# ============================================================================

datasets <- list(
    human   = "HumanRawCount_matched.csv",
    mouse   = "MouseRawCount_matched.csv",
    control = "HumanControlRawCount_matched.csv"
)

# Same exclusion list as the main pipeline -- these samples are technical
# outliers regardless of which time window is being analyzed (e.g. mouse's
# X.m_152 falls inside the first 2 hours and is still excluded here for the
# same reason it was excluded from the full-timecourse run).
outliers_to_exclude <- list(
    human   = c("NBplusNog_450", "NBplusNog_530", "NBplusNog_290", "NBplusNog_370"),
    mouse   = c("X.m_152", "X.m_468"),
    control = c("X.E8_410")
)

# The early-window cutoff. Samples with an extracted time strictly greater
# than this are dropped, in addition to (not instead of) the outlier
# exclusion above.
MAX_TIME_MINUTES <- 120

# ----------------------------------------------------------------------------
# Extracts the timepoint (minutes) from a column name. Identical logic to
# 01_diagnostics.R -- see that script for the full rationale/examples.
# ----------------------------------------------------------------------------
extract_time <- function(x) {
    m <- regmatches(x, regexpr("_[0-9]+", x))
    if (length(m) == 0 || m == "") return(NA_real_)
    as.numeric(sub("_", "", m))
}

# ============================================================================
# STAGE 1: Load the CSV, clean it, exclude outlier samples, AND restrict to
# samples with extracted time <= max_time. The time filter is applied last,
# after outlier exclusion, so an outlier sample that happens to fall inside
# the early window is still removed (not accidentally kept because it's
# "early").
# ============================================================================
load_and_clean <- function(csv_path, label, exclude_samples = character(0),
                            max_time = Inf) {
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

    if (is.finite(max_time)) {
        sample_time <- vapply(colnames(count_data), extract_time, numeric(1))
        keep_time <- !is.na(sample_time) & sample_time <= max_time
        cat(
            label, "restricting to <=", max_time, "minutes: keeping",
            sum(keep_time), "of", length(keep_time), "remaining sample(s)\n"
        )
        if (any(is.na(sample_time))) {
            warning(
                label, ": could not extract a time value for ",
                sum(is.na(sample_time)), " column(s) -- these were dropped ",
                "by the time filter. Check the naming pattern if unexpected."
            )
        }
        count_data <- count_data[, keep_time, drop = FALSE]
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
# Identical to 01_diagnostics.R. Normalization is recomputed from scratch
# on whatever count_matrix it's given, which here is already restricted to
# the early window.
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
# STAGE 3: Identical to 01_diagnostics.R. Writes CSV/PDF files under
# OUTPUT_ROOT/<dataset>/ (i.e. "1_Output_120min"/<dataset>/ here).
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
            main = paste(label, "(<=120min): VST-Normalized Expression by Sample (ordered by time)"),
            ylab = "VST value",
            col = time_colors,
            cex.axis = 0.6)
    abline(h = median(vst_ordered), col = "black", lty = 2)
    dev.off()

    pdf(out_path(label, "VST_Mean_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.25), height = 6)
    plot(time_ordered, sample_stats$mean,
         type = "b", pch = 19, col = time_colors,
         main = paste(label, "(<=120min): Mean VST Value Over Time"),
         xlab = "Time (minutes)", ylab = "Mean VST value")
    abline(h = mean(sample_stats$mean), col = "grey40", lty = 2)
    dev.off()

    sample_cor <- cor(vst_ordered, method = "pearson")
    write.csv(sample_cor, file = out_path(label, "VST_Sample_Correlation_Matrix.csv"))

    pdf(out_path(label, "VST_Sample_Correlation_Heatmap.pdf"), width = 7, height = 7)
    heatmap(sample_cor,
            symm = TRUE,
            Rowv = NA, Colv = NA,
            main = paste(label, "(<=120min): Sample-Sample Correlation (VST, ordered by time)"),
            margins = c(8, 8))
    dev.off()

    n_samples <- ncol(sample_cor)
    cor_list <- lapply(seq_len(n_samples), function(i) sample_cor[-i, i])
    names(cor_list) <- colnames(sample_cor)

    pdf(out_path(label, "VST_Correlation_Boxplot_By_Sample.pdf"),
        width = max(7, n_samples * 0.25), height = 6)
    boxplot(cor_list,
            las = 2,
            main = paste(label, "(<=120min): Per-Sample Correlation to All Other Samples"),
            ylab = "Correlation with other samples",
            col = time_colors,
            cex.axis = 0.6)
    abline(h = median(unlist(cor_list)), col = "black", lty = 2)
    dev.off()

    cat(label, "QC files written to ", OUTPUT_ROOT, "/", dataset_folder_names[[label]], "/\n--- End QC ---\n\n", sep = "")

    vst_ordered
}

# ============================================================================
# STAGE 2b: Gene-level mean expression histogram. Identical to
# 01_diagnostics.R. See that script for full rationale.
# ============================================================================
plot_mean_expression_histogram <- function(dds, label, sd_threshold = 3) {
    norm_counts <- counts(dds, normalized = TRUE)
    gene_means <- rowMeans(norm_counts)
    log_means <- log10(gene_means + 1)

    pdf(out_path(label, "Gene_Mean_Expression_Histogram.pdf"), width = 8, height = 6)
    hist(log_means,
         breaks = 60,
         col = "steelblue",
         border = "white",
         main = paste(label, "(<=120min): Per-Gene Mean Expression"),
         xlab = "log10(mean normalized count + 1)")
    abline(v = mean(log_means), col = "black", lty = 2)
    abline(v = mean(log_means) + sd_threshold * sd(log_means), col = "red", lty = 2)
    abline(v = mean(log_means) - sd_threshold * sd(log_means), col = "red", lty = 2)
    legend("topright", legend = c("mean", paste0("\u00b1", sd_threshold, " SD")),
           col = c("black", "red"), lty = 2, bty = "n")
    dev.off()
    cat(label, "wrote Gene_Mean_Expression_Histogram.pdf\n")

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
# STAGE 4: Gene filtering. Identical logic to 01_diagnostics.R (mean > 1,
# then top-quartile CV by this SUBSET's own distribution
# recomputed fresh here, since fewer samples generally means noisier CV estimates.
# It's important this isn't just reusing the full-timecourse thresholds).
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
    hist(cv_mean_passed,
         breaks = 60,
         col = "darkolivegreen",
         border = "white",
         main = paste(label, "(<=120min): Per-Gene CV (mean > 1 genes only)"),
         xlab = "Coefficient of variation (sd / mean, on VST scale)")
    abline(v = quantile(cv_mean_passed, 0.75), col = "red", lty = 2)
    legend("topright", legend = "75th percentile", col = "red", lty = 2, bty = "n")
    dev.off()
    cat(label, "wrote Gene_CV_Histogram.pdf\n")

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
# STAGE 5a: Cluster overlap evaluation
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
        title(main = paste(label, "(<=120min): Cluster Overlap for c =", k))
    }
    dev.off()

    m_opt
}

# ============================================================================
# ORCHESTRATOR
# ============================================================================
run_diagnostics <- function(csv_path, label, exclude_samples = character(0),
                             max_time = Inf) {
    cat("\n================ ", label, " (<=", max_time, " min) ================\n", sep = "")

    count_matrix <- load_and_clean(
        csv_path, label,
        exclude_samples = exclude_samples,
        max_time = max_time
    )

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
        exclude_samples = outliers_to_exclude[[label]],
        max_time = MAX_TIME_MINUTES
    )
}

cat("\nAll datasets diagnosed (<=", MAX_TIME_MINUTES, " min):", paste(names(datasets), collapse = ", "), "\n", sep = "")
cat("Files are under \"", OUTPUT_ROOT, "\"/<dataset>/. Next: review each dataset's\n", sep = "")
cat("Mfuzz_Cluster_Overlap_Evaluation.pdf, choose a c_final, then run 02_clustering_120min.R.\n")
