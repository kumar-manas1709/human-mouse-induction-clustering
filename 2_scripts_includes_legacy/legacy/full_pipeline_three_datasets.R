library(readr)
library(DESeq2)
library(Mfuzz)
library(Biobase)

datasets <- list(
    human   = "HumanRawCount_matched.csv",
    mouse   = "MouseRawCount_matched.csv",
    control = "HumanControlRawCount_matched.csv"
)

outliers_to_exclude <- list(
    human   = c("NBplusNog_450", "NBplusNog_530"),
    mouse   = character(0),
    control = c("X.E8_410")
)

c_final_choices <- list(
    human   = 6,
    mouse   = 6,
    control = 6
)
extract_time <- function(x) {
    m <- regmatches(x, regexpr("_[0-9]+", x))
    if (length(m) == 0 || m == "") return(NA_real_)
    as.numeric(sub("_", "", m))
}

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
        sample = colnames(vst_ordered),
        time_minutes = time_ordered,
        mean   = colMeans(vst_ordered),
        median = apply(vst_ordered, 2, median),
        sd     = apply(vst_ordered, 2, sd)
    )
    print(sample_stats)
    write_csv(sample_stats, paste0(label, "_VST_Sample_Summary_Stats.csv"))

    pdf(paste0(label, "_VST_Boxplot_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.25), height = 6)
    boxplot(vst_ordered,
            las = 2,
            main = paste(label, ": VST-Normalized Expression by Sample (ordered by time)"),
            ylab = "VST value",
            col = time_colors,
            cex.axis = 0.6)
    abline(h = median(vst_ordered), col = "black", lty = 2)
    dev.off()

    pdf(paste0(label, "_VST_Mean_By_Sample.pdf"), width = max(7, ncol(vst_ordered) * 0.25), height = 6)
    plot(time_ordered, sample_stats$mean,
         type = "b", pch = 19, col = time_colors,
         main = paste(label, ": Mean VST Value Over Time"),
         xlab = "Time (minutes)", ylab = "Mean VST value")
    abline(h = mean(sample_stats$mean), col = "grey40", lty = 2)
    dev.off()

    sample_cor <- cor(vst_ordered, method = "pearson")
    write.csv(sample_cor, file = paste0(label, "_VST_Sample_Correlation_Matrix.csv"))

    pdf(paste0(label, "_VST_Sample_Correlation_Heatmap.pdf"), width = 7, height = 7)
    heatmap(sample_cor,
            symm = TRUE,
            Rowv = NA, Colv = NA,
            main = paste(label, ": Sample-Sample Correlation (VST, ordered by time)"),
            margins = c(8, 8))
    dev.off()

    cat(label, "QC files written.\n--- End QC ---\n\n")

    vst_ordered
}
filter_genes <- function(dds, vst_matrix, label) {
    norm_counts <- counts(dds, normalized = TRUE)
    keep_mean <- rowMeans(norm_counts) > 1
    cat(label, "genes with mean > 1:", sum(keep_mean), "\n")

    cv <- apply(vst_matrix, 1, sd) / abs(rowMeans(vst_matrix))
    cv_threshold <- 0.4
    keep_cv <- cv > cv_threshold
    cat(label, "genes with CV >", cv_threshold, ":", sum(keep_cv), "\n")

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

evaluate_cluster_overlap <- function(eset_r, label) {
    m_opt <- mestimate(eset_r)
    cat(label, "m_opt:", m_opt, "\n")
    pdf(paste0(label, "_Mfuzz_Cluster_Overlap_Evaluation.pdf"), width = 7, height = 7)
    for (k in c(4, 6, 8, 10)) {
        set.seed(42)
        cl_test <- mfuzz(eset_r, c = k, m = m_opt)
        O <- overlap(cl_test)
        overlap.plot(cl_test, over = O, thres = 0.1)
        title(main = paste(label, ": Cluster Overlap for c =", k))
    }
    dev.off()

    m_opt
}

preview_cluster_counts <- function(eset_r, m_opt, label, candidates = c(4, 6, 8, 10)) {
    for (k in candidates) {
        set.seed(42)
        cl_preview <- mfuzz(eset_r, c = k, m = m_opt)
        preview_label <- paste0(label, "_c", k, "_preview")
        plot_cluster_centroids(cl_preview, preview_label)
    }
    cat(label, ": wrote centroid previews for c =", paste(candidates, collapse = ", "), "\n")
}

plot_cluster_centroids <- function(final_clusters, label) {
    centers <- final_clusters$centers   # rows = clusters, cols = samples
    c_final <- nrow(centers)

    okabe_ito <- c(
        "#E69F00", # orange
        "#56B4E9", # sky blue
        "#009E73", # bluish green
        "#F0E442", # yellow
        "#0072B2", # blue
        "#D55E00", # vermillion
        "#CC79A7", # reddish purple
        "#000000"  # black
    )
    cluster_colors <- rep_len(okabe_ito, c_final)

    cluster_ltys <- rep_len(1:6, c_final)

    png(paste0(label, "_Mfuzz_Cluster_Centroids.png"), width = 1400, height = 900, res = 150)
    par(mar = c(8, 4, 4, 8), xpd = TRUE)  # extra margin for rotated labels + legend
    matplot(t(centers),
            type = "l", lty = cluster_ltys, lwd = 2.5, col = cluster_colors,
            xaxt = "n",
            main = paste(label, ": Cluster Centroids (mean trajectory per cluster)"),
            xlab = "", ylab = "Mean standardized expression")
    axis(1, at = seq_len(ncol(centers)), labels = colnames(centers), las = 2, cex.axis = 0.5)
    mtext("Time", side = 1, line = 6.5)
    legend("topright", inset = c(-0.15, 0),
           legend = paste("Cluster", seq_len(c_final)),
           col = cluster_colors, lty = cluster_ltys, lwd = 2.5, bty = "n", cex = 0.8)
    dev.off()
    cat(label, "wrote", paste0(label, "_Mfuzz_Cluster_Centroids.png"), "\n")
}
run_final_clustering <- function(eset_r, final_vst_data, m_opt, c_final, label) {
    set.seed(42)
    final_clusters <- mfuzz(eset_r, c = c_final, m = m_opt)

    grid_rows <- ceiling(sqrt(c_final))
    grid_cols <- ceiling(c_final / grid_rows)

    png(paste0(label, "_Mfuzz_Final_Clusters_Plot.png"),
        width = 1200, height = 1000, res = 150)
    mfuzz.plot2(
        eset_r,
        cl = final_clusters,
        mfrow = c(grid_rows, grid_cols),
        time.labels = colnames(exprs(eset_r)),
        x11 = FALSE
    )
    dev.off()
    cat(label, "wrote", paste0(label, "_Mfuzz_Final_Clusters_Plot.png"), "\n")

    plot_cluster_centroids(final_clusters, label)

    gene_assignments <- final_clusters$cluster
    membership_matrix <- final_clusters$membership

    clustering_results <- data.frame(
        Gene = rownames(final_vst_data),
        Assigned_Cluster = gene_assignments,
        membership_matrix
    )
    clustering_results <- clustering_results[order(clustering_results$Assigned_Cluster), ]
    write_csv(clustering_results, paste0(label, "_Mfuzz_Soft_Clustering_Results.csv"))

    list(clusters = final_clusters, clustering_results = clustering_results)
}

run_pipeline <- function(csv_path, label, c_final = 6, exclude_samples = character(0),
                          preview_counts = TRUE) {
    cat("\n================ ", label, " ================\n")

    count_matrix <- load_and_clean(csv_path, label, exclude_samples = exclude_samples)

    norm <- normalize_deseq2(count_matrix)
    dds <- norm$dds
    vst_matrix <- norm$vst_matrix

    run_vst_qc(vst_matrix, label)

    final_vst_data <- filter_genes(dds, vst_matrix, label)

    eset <- ExpressionSet(assayData = as.matrix(final_vst_data))
    eset_r <- standardise(eset)

    m_opt <- evaluate_cluster_overlap(eset_r, label)

    if (preview_counts) {
        preview_cluster_counts(eset_r, m_opt, label)
    }

    cluster_out <- run_final_clustering(eset_r, final_vst_data, m_opt, c_final, label)

    cat(label, ": pipeline complete.\n")

    invisible(list(
        vst_matrix = vst_matrix,
        final_vst_data = final_vst_data,
        clusters = cluster_out$clusters,
        clustering_results = cluster_out$clustering_results
    ))
}
results <- list()
for (label in names(datasets)) {
    results[[label]] <- run_pipeline(
        datasets[[label]], label, c_final = c_final_choices[[label]],
        exclude_samples = outliers_to_exclude[[label]]
    )
}

cat("\nAll datasets processed:", paste(names(datasets), collapse = ", "), "\n")
