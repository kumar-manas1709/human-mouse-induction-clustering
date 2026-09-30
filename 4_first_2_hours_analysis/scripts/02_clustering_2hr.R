library(readr)
library(Mfuzz)
library(Biobase)

# ============================================================================
# 02_clustering_2hr.R
#
# Variant of ../../scripts/02_clustering.R for the <=120-minute subset.
# Picks up from 1_Output_First2Hours/<dataset>/pipeline_state.rds written
# by 01_diagnostics_2hr.R -- run that first.
#
# USAGE
# -----
#   Rscript 02_clustering_2hr.R
#
# (run from this folder -- first_2_hours_analysis/scripts/)
#
# OUTPUT
# ------
# Previews, per dataset, under 1_Output_First2Hours/<dataset>/:
#   Mfuzz_Cluster_Centroids_c{2,3,4,6,8,10}_preview.png
# Final deliverables, shared folder, under 1_Output_First2Hours/4_clusters/:
#   <dataset>_Mfuzz_Final_Clusters_Plot.png
#   <dataset>_Mfuzz_Cluster_Centroids.png
#   <dataset>_Mfuzz_Soft_Clustering_Results.csv
# ============================================================================

MAX_TIME_MINUTES <- 120

labels <- c("human", "mouse", "control")

# Set to FALSE to skip the interactive pause and just use c_final_choices.
prompt_for_c_final <- TRUE
c_final_choices <- list(human = 3, mouse = 3, control = 3)
preview_counts <- TRUE

dataset_folder_names <- c(human = "1_human", mouse = "2_mouse", control = "3_control")

out_path <- function(dataset, filename) {
    folder_name <- if (dataset %in% names(dataset_folder_names)) dataset_folder_names[[dataset]] else dataset
    dir <- file.path("1_Output_First2Hours", folder_name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file.path(dir, filename)
}

final_output_path <- function(dataset, filename) {
    dir <- file.path("1_Output_First2Hours", "4_clusters")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file.path(dir, paste0(dataset, "_", filename))
}

preview_cluster_counts <- function(eset_r, m_opt, label, candidates = c(2, 3, 4, 6, 8, 10)) {
    for (k in candidates) {
        set.seed(42)
        cl_preview <- mfuzz(eset_r, c = k, m = m_opt)
        plot_cluster_centroids(
            cl_preview, save_path = out_path(label, paste0("Mfuzz_Cluster_Centroids_c", k, "_preview.png")),
            title_label = paste0(label, " (c=", k, ", <=", MAX_TIME_MINUTES, "min, preview)")
        )
    }
    cat(label, ": wrote centroid previews for c =", paste(candidates, collapse = ", "), "\n")
}

plot_cluster_centroids <- function(final_clusters, save_path, title_label) {
    centers <- final_clusters$centers
    c_final <- nrow(centers)

    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    cluster_colors <- rep_len(okabe_ito, c_final)
    cluster_ltys <- rep_len(1:6, c_final)

    png(save_path, width = 1400, height = 900, res = 150)
    par(mar = c(8, 4, 4, 8), xpd = TRUE)
    matplot(t(centers), type = "l", lty = cluster_ltys, lwd = 2.5, col = cluster_colors, xaxt = "n",
            main = paste(title_label, ": Cluster Centroids"), xlab = "", ylab = "Mean standardized expression")
    axis(1, at = seq_len(ncol(centers)), labels = colnames(centers), las = 2, cex.axis = 0.6)
    mtext("Time", side = 1, line = 6.5)
    legend("topright", inset = c(-0.15, 0), legend = paste("Cluster", seq_len(c_final)),
           col = cluster_colors, lty = cluster_ltys, lwd = 2.5, bty = "n", cex = 0.8)
    dev.off()
    cat("wrote", save_path, "\n")
}

run_final_clustering <- function(eset_r, final_vst_data, m_opt, c_final, label) {
    set.seed(42)
    final_clusters <- mfuzz(eset_r, c = c_final, m = m_opt)

    grid_rows <- ceiling(sqrt(c_final))
    grid_cols <- ceiling(c_final / grid_rows)
    final_plot_path <- final_output_path(label, "Mfuzz_Final_Clusters_Plot.png")
    png(final_plot_path, width = 1200, height = 1000, res = 150)
    mfuzz.plot2(eset_r, cl = final_clusters, mfrow = c(grid_rows, grid_cols),
                time.labels = colnames(exprs(eset_r)), x11 = FALSE)
    dev.off()
    cat("wrote", final_plot_path, "\n")

    plot_cluster_centroids(final_clusters, save_path = final_output_path(label, "Mfuzz_Cluster_Centroids.png"), title_label = label)

    clustering_results <- data.frame(
        Dataset = label, Gene = rownames(final_vst_data),
        Assigned_Cluster = final_clusters$cluster, final_clusters$membership
    )
    clustering_results <- clustering_results[order(clustering_results$Assigned_Cluster), ]
    csv_path <- final_output_path(label, "Mfuzz_Soft_Clustering_Results.csv")
    write_csv(clustering_results, csv_path)
    cat("wrote", csv_path, "\n")

    list(clusters = final_clusters, clustering_results = clustering_results)
}

get_c_final_interactively <- function(label, default) {
    prompt <- paste0(
        "\n", label, " (<=", MAX_TIME_MINUTES, "min): review 1_Output_First2Hours/", label,
        "/Mfuzz_Cluster_Overlap_Evaluation.pdf and the c*_preview centroid files just written.\n",
        "Enter number of clusters to use [default ", default, "]: "
    )
    repeat {
        if (interactive()) input <- readline(prompt) else { cat(prompt); input <- readLines("stdin", n = 1) }
        if (identical(input, "")) { cat(label, ": using default c_final =", default, "\n"); return(default) }
        val <- suppressWarnings(as.integer(input))
        if (!is.na(val) && val > 1) return(val)
        cat("Please enter a whole number greater than 1 (or press Enter for the default).\n")
    }
}

run_clustering <- function(label, default_c_final, preview_counts = TRUE, prompt_for_c_final = TRUE) {
    state_path <- out_path(label, "pipeline_state.rds")
    if (!file.exists(state_path)) {
        stop("Could not find ", state_path, " -- run 01_diagnostics_2hr.R first.")
    }
    state <- readRDS(state_path)
    if (preview_counts) preview_cluster_counts(state$eset_r, state$m_opt, label)
    c_final <- if (prompt_for_c_final) get_c_final_interactively(label, default_c_final) else default_c_final
    cluster_out <- run_final_clustering(state$eset_r, state$final_vst_data, state$m_opt, c_final, label)
    cat(label, ": clustering complete (c_final =", c_final, ").\n")
    invisible(list(final_vst_data = state$final_vst_data, clusters = cluster_out$clusters,
                    clustering_results = cluster_out$clustering_results))
}

results <- list()
for (label in labels) {
    results[[label]] <- run_clustering(label, default_c_final = c_final_choices[[label]],
                                        preview_counts = preview_counts, prompt_for_c_final = prompt_for_c_final)
}

cat("\nAll datasets clustered (<=", MAX_TIME_MINUTES, "min):", paste(labels, collapse = ", "), "\n")
