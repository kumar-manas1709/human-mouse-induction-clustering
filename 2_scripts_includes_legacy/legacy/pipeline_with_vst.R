library(readODS)
library(DESeq2)
library(Mfuzz)
library(Biobase)

count_data <- read_ods("RawData/HumanRawNew.ods", sheet = 1, as_tibble = FALSE, row_names = TRUE)
dim(count_data) # Dimension verification

count_data[is.na(count_data)] <- 0
valid_rows <- !is.na(rownames(count_data)) & rownames(count_data) != ""
count_data <- count_data[valid_rows, ]
count_matrix <- as.matrix(count_data)
count_matrix <- round(count_matrix)

count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
nrow(count_matrix) # Number of genes left

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

extract_last_number <- function(x) {
    matches <- regmatches(x, gregexpr("[0-9]+", x))[[1]]
    if (length(matches) == 0) return(NA_real_)
    as.numeric(matches[length(matches)])
}
sample_time <- vapply(colnames(vst_matrix), extract_last_number, numeric(1))
if (any(is.na(sample_time))) {
    warning(
        "Could not extract a trailing time value from these column names: ",
        paste(colnames(vst_matrix)[is.na(sample_time)], collapse = ", "),
        " -- check the naming pattern. Falling back to original column order."
    )
    sample_time <- seq_len(ncol(vst_matrix))
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
write.csv(sample_stats, file = "VST_Sample_Summary_Stats.csv", row.names = FALSE)
pdf("VST_Boxplot_By_Sample.pdf", width = max(7, ncol(vst_ordered) * 0.4), height = 6)
boxplot(vst_ordered,
        las = 2,                     # rotate sample labels for readability
        main = "VST-Normalized Expression by Sample (ordered by time)",
        ylab = "VST value",
        col = time_colors,
        cex.axis = 0.7)
abline(h = median(vst_ordered), col = "black", lty = 2)  # overall median for reference
dev.off()


pdf("VST_Mean_By_Sample.pdf", width = max(7, ncol(vst_ordered) * 0.4), height = 6)
plot(time_ordered, sample_stats$mean,
     type = "b", pch = 19, col = time_colors,
     main = "Mean VST Value Over Time",
     xlab = "Time (minutes)", ylab = "Mean VST value")
abline(h = mean(sample_stats$mean), col = "grey40", lty = 2)
dev.off()

sample_cor <- cor(vst_ordered, method = "pearson")
write.csv(sample_cor, file = "VST_Sample_Correlation_Matrix.csv")

pdf("VST_Sample_Correlation_Heatmap.pdf", width = 7, height = 7)
heatmap(sample_cor,
        symm = TRUE,
        Rowv = NA, Colv = NA,   # keep time order instead of re-clustering
        main = "Sample-Sample Correlation (VST, ordered by time)",
        margins = c(8, 8))
dev.off()

norm_counts <- counts(dds, normalized = TRUE)
keep_mean <- rowMeans(norm_counts) > 1
cat("Genes with mean > 1:", sum(keep_mean), "\n")

cv <- apply(vst_matrix, 1, sd) / abs(rowMeans(vst_matrix))
cv_threshold <- 0.4
keep_cv <- cv > cv_threshold
cat("Genes with CV >", cv_threshold, ":", sum(keep_cv), "\n")

final_vst_data <- vst_matrix[keep_mean & keep_cv, ]
cat("Genes passing BOTH filters:", nrow(final_vst_data), "\n")

if (nrow(final_vst_data) == 0) {
    warning("CV filter removed all genes! Falling back to top 2,000 variable genes (with mean > 1).")
    vst_high_mean <- vst_matrix[keep_mean, ]
    gene_var <- apply(vst_high_mean, 1, var)
    top_n <- min(2000, nrow(vst_high_mean)) 
    keep_top <- order(gene_var, decreasing = TRUE)[1:top_n]
    final_vst_data <- vst_high_mean[keep_top, ]
    cat("Fallback: retained", nrow(final_vst_data), "genes by variance.\n")
}

# Incase you want to use the filtered data elsewhere
# cat("EXPORTING FILTERED VST TO CSV...\n")
# write.csv(as.data.frame(final_vst_data), 
#          file = "vst_filtered_for_clustering.csv", 
#          row.names = TRUE)

eset <- ExpressionSet(assayData = as.matrix(final_vst_data))
eset_r <- standardise(eset)


m_opt <- mestimate(eset_r)
m_opt

# The overlap is done twice for my own understanding. Not needed, only the more sensitive one is necessary!
pdf("Mfuzz_Cluster_Overlap_Evaluation.pdf", width = 7, height = 7)
for (k in c(4, 6, 8, 10)) {
    cl_test <- mfuzz(eset_r, c = k, m = m_opt)
    O <- overlap(cl_test)
    overlap.plot(cl_test, over = O, thres = 0.5)
    title(main = paste("Cluster Overlap for c =", k)) 
            
}
dev.off()

pdf("Mfuzz_Cluster_Overlap_Evaluation_Stable.pdf", width = 7, height = 7)
for (k in c(4, 6, 8, 10)) {
    set.seed(42) 
    cl_test <- mfuzz(eset_r, c = k, m = m_opt)
    O <- overlap(cl_test)
    overlap.plot(cl_test, over = O, thres = 0.1)
    title(main = paste("Cluster Overlap for c =", k))
}
dev.off()

# DO NOT RUN THE FOLLOWING CODE IN INTERACTIVE!!! ONLY IN SOURCE!!!
# {
#   cluster_options <- c(4, 6, 8, 10)
#   menu_labels <- paste(cluster_options, "clusters")
#   
#   cat("\n--- FINAL CLUSTER SELECTION ---\n")
#   user_selection <- menu(menu_labels, title = "Select your final number of clusters:")
#  
#   if (user_selection == 0) {
#       stop("Pipeline execution cancelled.")
#   } else {
#       c_final <- cluster_options[user_selection]
#   }
# }
# IN INTERACTIVE, STOP HERE AND THEN MANUALLY ASSIGN THE CLUSTER NUMBER!!! PLEASE!!!

c_final <- 6

set.seed(42) 
final_clusters <- mfuzz(eset_r, c = c_final, m = m_opt)
set.seed(42) 
final_clusters <- mfuzz(eset_r, c = c_final, m = m_opt)

grid_rows <- ceiling(sqrt(c_final))
grid_cols <- ceiling(c_final / grid_rows)
png("Mfuzz_Final_Clusters_Plot.png", width = 1200, height = 1000, res = 150)
dev.off()

gene_assignments <- final_clusters$cluster
membership_matrix <- final_clusters$membership

clustering_results <- data.frame(
    Gene = rownames(final_vst_data),
    Assigned_Cluster = gene_assignments,
    membership_matrix
)

clustering_results <- clustering_results[order(clustering_results$Assigned_Cluster), ]
write.csv(clustering_results, 
          file = "Mfuzz_Soft_Clustering_Results.csv", 
          row.names = FALSE)
