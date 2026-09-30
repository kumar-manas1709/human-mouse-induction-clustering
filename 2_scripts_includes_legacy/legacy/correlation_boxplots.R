
suppressMessages(library(readr))

source("pipeline_utils.R")  # provides out_path(dataset, filename)

extract_time <- function(x) {
    m <- regmatches(x, regexpr("_[0-9]+", x))
    if (length(m) == 0 || m == "") return(NA_real_)
    as.numeric(sub("_", "", m))
}

plot_correlation_boxplot <- function(label) {
    cat("Processing", label, "...\n")

    matrix_path <- out_path(label, "VST_Sample_Correlation_Matrix.csv")
    if (!file.exists(matrix_path)) {
        warning(label, ": could not find ", matrix_path, " -- run 01_diagnostics.R first. Skipping.")
        return(invisible(NULL))
    }

    cor_mat <- read.csv(matrix_path, row.names = 1, check.names = FALSE)
    cor_mat <- as.matrix(cor_mat)

    # Order samples by time, same convention as the rest of the pipeline.
    sample_time <- vapply(colnames(cor_mat), extract_time, numeric(1))
    time_order <- order(sample_time)
    cor_mat <- cor_mat[time_order, time_order]
    time_ordered <- sample_time[time_order]

    # For each sample (column), take its correlation with every OTHER sample
    # -- i.e. the column with the diagonal (self=1) removed.
    n <- ncol(cor_mat)
    cor_list <- lapply(seq_len(n), function(i) cor_mat[-i, i])
    names(cor_list) <- colnames(cor_mat)

    time_colors <- colorRampPalette(c("blue", "red"))(length(unique(time_ordered)))[
        as.integer(factor(time_ordered))
    ]

    pdf(out_path(label, "VST_Correlation_Boxplot_By_Sample.pdf"),
        width = max(7, n * 0.25), height = 6)
    boxplot(cor_list,
            las = 2,
            main = paste(label, ": Per-Sample Correlation to All Other Samples"),
            ylab = "Correlation with other samples",
            col = time_colors,
            cex.axis = 0.6)
    abline(h = median(unlist(cor_list)), col = "black", lty = 2)
    dev.off()

    cat("  wrote VST_Correlation_Boxplot_By_Sample.pdf to Output/", label, "/\n", sep = "")
}

labels <- c("human", "mouse", "control")

for (label in labels) {
    plot_correlation_boxplot(label)
}

cat("\nDone.\n")
