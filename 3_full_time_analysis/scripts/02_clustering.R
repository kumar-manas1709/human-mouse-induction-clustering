library(readr)
library(Mfuzz)
library(Biobase)

source("pipeline_utils.R")  # provides out_path(dataset, filename) and final_output_path(dataset, filename)

# ============================================================================
# 02_clustering.R
#
# Part 2 of 2: centroid previews across candidate cluster counts, then a
# PAUSE where you're asked to type in the number of clusters to use, then
# final Mfuzz clustering.
#
# Picks up from the "1_Output"/<dataset>/pipeline_state.rds files written by
# 01_diagnostics.R -- run that script first. Because this script reads
# saved state instead of re-normalizing, you can re-run it as many times as
# you like with different cluster counts without waiting on DESeq2/VST
# again.
#
# USAGE
# -----
#   Rscript 02_clustering.R
#
# (run from the same folder as 01_diagnostics.R was run from, so it can
# find pipeline_utils.R and the "1_Output"/<dataset>/ folders it created)
#
# For each dataset in turn, this will:
#   1. write the centroid preview PNGs for c = 2, 3, 4, 6, 8, 10
#   2. PAUSE and print a prompt in the console asking you to enter a
#      cluster count for that dataset (go look at the files from step 1,
#      then type a number and press Enter -- or just press Enter to accept
#      the suggested default)
#   3. run final clustering with whatever number you entered
# ...before moving on to the next dataset. Set prompt_for_c_final = FALSE
# below to skip the pause entirely and just use c_final_choices instead.
#
# OUTPUT
# ------
# The EXPLORATORY previews (still per-dataset, since you're comparing
# candidate c values one dataset at a time) go under
# "1_Output"/<numbered dataset folder>/:
#   Mfuzz_Cluster_Centroids_c2_preview.png
#   Mfuzz_Cluster_Centroids_c3_preview.png
#   Mfuzz_Cluster_Centroids_c4_preview.png
#   Mfuzz_Cluster_Centroids_c6_preview.png
#   Mfuzz_Cluster_Centroids_c8_preview.png
#   Mfuzz_Cluster_Centroids_c10_preview.png
#
# The FINAL deliverables -- once you've committed to a c_final -- go into
# one shared folder, "1_Output"/"4_clusters"/, with the dataset name
# prepended to each filename so all three datasets' finished results sit
# together for easy side-by-side comparison:
#   human_Mfuzz_Final_Clusters_Plot.png     (and mouse_..., control_...)
#   human_Mfuzz_Cluster_Centroids.png       (and mouse_..., control_...)
#   human_Mfuzz_Soft_Clustering_Results.csv (and mouse_..., control_...)
# ============================================================================

# ============================================================================
# CONFIG
# ============================================================================

labels <- c("human", "mouse", "control")

# If TRUE, the script pauses after writing each dataset's overlap PDF and
# centroid previews, and asks you to type in the number of clusters to use
# for that dataset -- go look at the files it just wrote, then answer in
# the console/terminal. If FALSE, it uses c_final_choices below instead
# (no pause) -- handy once you already know your final numbers and just
# want to (re-)run everything in one go, e.g. non-interactively.
prompt_for_c_final <- TRUE

# Only used when prompt_for_c_final = FALSE, and as the suggested default
# shown in the prompt when prompt_for_c_final = TRUE (press Enter to accept
# it without typing a number).
c_final_choices <- list(
    human   = 6,
    mouse   = 6,
    control = 6
)

# Set to FALSE once you've settled on a c_final for each dataset -- these
# previews are only useful while you're still deciding, and re-running
# Mfuzz for every candidate c adds runtime.
preview_counts <- TRUE

# ============================================================================
# STAGE 5a-helper: For each candidate c value, run Mfuzz and immediately
# plot its centroids -- gives you a direct, visual side-by-side of cluster
# SHAPE at each candidate count, rather than only the overlap-matrix numbers.
# Look through these files for the smallest c where every cluster still
# looks visually distinct (no two lines tracking each other closely).
# ============================================================================
preview_cluster_counts <- function(eset_r, m_opt, label, candidates = c(2, 3, 4, 6, 8, 10)) {
    for (k in candidates) {
        set.seed(42)
        cl_preview <- mfuzz(eset_r, c = k, m = m_opt)
        plot_cluster_centroids(
            cl_preview,
            save_path = out_path(label, paste0("Mfuzz_Cluster_Centroids_c", k, "_preview.png")),
            title_label = paste0(label, " (c=", k, ", preview)")
        )
    }
    cat(label, ": wrote centroid previews for c =", paste(candidates, collapse = ", "), "\n")
}

# ============================================================================
# STAGE 5b-helper: Plot each cluster's CENTROID (mean trajectory) as one
# clean line per cluster, overlaid on a single plot. This is the readable
# alternative to mfuzz.plot2 when most genes have high membership and heavy
# overplotting turns the per-cluster panels into a solid block of color --
# it shows you the actual SHAPE each cluster represents, at a glance.
#
# Takes an explicit `save_path` rather than building one itself, so the
# caller decides WHERE it lands -- previews go under the per-dataset
# diagnostics folder (out_path), while the final centroid plot goes into
# the shared "4_clusters" folder (final_output_path) alongside the other
# final deliverables.
# ============================================================================
plot_cluster_centroids <- function(final_clusters, save_path, title_label) {
    centers <- final_clusters$centers   # rows = clusters, cols = samples
    c_final <- nrow(centers)

    # Okabe-Ito palette: designed to remain distinguishable under the most
    # common forms of color blindness (deuteranopia/protanopia/tritanopia),
    # unlike rainbow() where red and green (the two hardest colors to tell
    # apart for red-green color blindness) sit right next to each other.
    # Cycled with rep_len in case c_final > 8.
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

    # Distinct line types too, so clusters stay distinguishable even in
    # black-and-white printouts or for any residual color-perception
    # difficulty -- redundant coding (color + shape) rather than color alone.
    cluster_ltys <- rep_len(1:6, c_final)

    png(save_path, width = 1400, height = 900, res = 150)
    par(mar = c(8, 4, 4, 8), xpd = TRUE)  # extra margin for rotated labels + legend
    matplot(t(centers),
            type = "l", lty = cluster_ltys, lwd = 2.5, col = cluster_colors,
            xaxt = "n",
            main = paste(title_label, ": Cluster Centroids (mean trajectory per cluster)"),
            xlab = "", ylab = "Mean standardized expression")
    axis(1, at = seq_len(ncol(centers)), labels = colnames(centers), las = 2, cex.axis = 0.5)
    mtext("Time", side = 1, line = 6.5)
    legend("topright", inset = c(-0.15, 0),
           legend = paste("Cluster", seq_len(c_final)),
           col = cluster_colors, lty = cluster_ltys, lwd = 2.5, bty = "n", cex = 0.8)
    dev.off()
    cat("wrote", save_path, "\n")
}

# ============================================================================
# STAGE 5b: Run the final Mfuzz clustering at a chosen c_final, plot the
# cluster expression profiles (both the full per-gene panels AND the clean
# centroid overlay), and write the soft-clustering results CSV.
# ============================================================================
run_final_clustering <- function(eset_r, final_vst_data, m_opt, c_final, label) {
    set.seed(42)
    final_clusters <- mfuzz(eset_r, c = c_final, m = m_opt)

    # Cluster profile plot: one panel per cluster, showing every gene's
    # (standardized) expression trajectory across your samples, colored by
    # how strongly it belongs to that cluster (its membership score). With
    # many genes and high membership scores, this can look like a dense
    # block of color -- see the centroid plot below for a cleaner view of
    # cluster shape.
    grid_rows <- ceiling(sqrt(c_final))
    grid_cols <- ceiling(c_final / grid_rows)

    final_plot_path <- final_output_path(label, "Mfuzz_Final_Clusters_Plot.png")
    png(final_plot_path, width = 1200, height = 1000, res = 150)
    mfuzz.plot2(
        eset_r,
        cl = final_clusters,
        mfrow = c(grid_rows, grid_cols),
        time.labels = colnames(exprs(eset_r)),
        x11 = FALSE
    )
    dev.off()
    cat("wrote", final_plot_path, "\n")

    plot_cluster_centroids(
        final_clusters,
        save_path = final_output_path(label, "Mfuzz_Cluster_Centroids.png"),
        title_label = label
    )

    gene_assignments <- final_clusters$cluster
    membership_matrix <- final_clusters$membership

    clustering_results <- data.frame(
        Dataset = label,
        Gene = rownames(final_vst_data),
        Assigned_Cluster = gene_assignments,
        membership_matrix
    )
    clustering_results <- clustering_results[order(clustering_results$Assigned_Cluster), ]
    csv_path <- final_output_path(label, "Mfuzz_Soft_Clustering_Results.csv")
    write_csv(clustering_results, csv_path)
    cat("wrote", csv_path, "\n")

    list(clusters = final_clusters, clustering_results = clustering_results)
}

# ============================================================================
# Prompts for a cluster count in the console/terminal and waits for input --
# works whether this script is run via `Rscript` (reads from stdin) or
# interactively in RStudio/an R session (uses readline() instead). Keeps
# asking until you enter a whole number > 1, or just press Enter to accept
# the suggested default.
# ============================================================================
get_c_final_interactively <- function(label, default) {
    prompt <- paste0(
        "\n", label, ": review \"1_Output\"/", dataset_folder_names[[label]],
        "/Mfuzz_Cluster_Overlap_Evaluation.pdf ",
        "and the \"1_Output\"/", dataset_folder_names[[label]],
        "/Mfuzz_Cluster_Centroids_c*_preview.png files just written.\n",
        "Enter number of clusters to use for ", label, " [default ", default, "]: "
    )
    repeat {
        if (interactive()) {
            input <- readline(prompt)
        } else {
            cat(prompt)
            input <- readLines("stdin", n = 1)
        }
        if (identical(input, "")) {
            cat(label, ": using default c_final =", default, "\n")
            return(default)
        }
        val <- suppressWarnings(as.integer(input))
        if (!is.na(val) && val > 1) {
            return(val)
        }
        cat("Please enter a whole number greater than 1 (or press Enter for the default).\n")
    }
}

# ============================================================================
# ORCHESTRATOR: loads a dataset's saved state from 01_diagnostics.R, runs
# the centroid preview (optional), pauses for you to choose c_final (if
# prompt_for_c_final = TRUE), then runs final clustering.
# ============================================================================
run_clustering <- function(label, default_c_final, preview_counts = TRUE,
                            prompt_for_c_final = TRUE) {
    state_path <- out_path(label, "pipeline_state.rds")
    if (!file.exists(state_path)) {
        stop(
            "Could not find ", state_path, " -- run 01_diagnostics.R first ",
            "(it must be run from the same folder as this script)."
        )
    }
    state <- readRDS(state_path)

    if (preview_counts) {
        preview_cluster_counts(state$eset_r, state$m_opt, label)
    }

    c_final <- if (prompt_for_c_final) {
        get_c_final_interactively(label, default_c_final)
    } else {
        default_c_final
    }

    cluster_out <- run_final_clustering(
        state$eset_r, state$final_vst_data, state$m_opt, c_final, label
    )

    cat(label, ": clustering complete (c_final =", c_final, ").\n")

    invisible(list(
        final_vst_data = state$final_vst_data,
        clusters = cluster_out$clusters,
        clustering_results = cluster_out$clustering_results
    ))
}

# ============================================================================
# Run clustering once per dataset
# ============================================================================

results <- list()
for (label in labels) {
    results[[label]] <- run_clustering(
        label,
        default_c_final = c_final_choices[[label]],
        preview_counts = preview_counts,
        prompt_for_c_final = prompt_for_c_final
    )
}

cat("\nAll datasets clustered:", paste(labels, collapse = ", "), "\n")
cat("Files are under \"1_Output\"/<dataset>/.\n")
