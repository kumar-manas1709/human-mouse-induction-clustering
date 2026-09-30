# ============================================================================
# pipeline_utils.R
#
# Shared helper sourced by both 01_diagnostics.R and 02_clustering.R.
# Keep this file in the same folder as both scripts.
# ============================================================================

# Numbered folder names, purely so the three dataset folders sort in this
# specific order (human, mouse, control) in a file browser instead of
# whatever order alphabetical sorting would give ("control" would otherwise
# sort before "human"). Add an entry here if a new dataset label is ever
# introduced. Any label not listed falls back to its plain name with no
# number prefix.

dataset_folder_names <- c(
    human   = "1_human",
    mouse   = "2_mouse",
    control = "3_control"
)

# Base output folder. Defaults to "1_Output" for the normal full-timecourse
# pipeline. A script can point its outputs somewhere else entirely (e.g. a
# parallel analysis on a subset of samples) by setting OUTPUT_ROOT to a
# different value BEFORE calling source("pipeline_utils.R")

if (!exists("OUTPUT_ROOT")) OUTPUT_ROOT <- "1_Output"

# Builds a path like "1_Output/1_human/VST_Boxplot_By_Sample.pdf",
# creating the folder first if it doesn't already exist. Every output file
# from both scripts is written through this, so everything ends up grouped
# by dataset instead of all mixed into one flat folder.
out_path <- function(dataset, filename) {
    folder_name <- if (dataset %in% names(dataset_folder_names)) {
        dataset_folder_names[[dataset]]
    } else {
        dataset  # fallback for any label not in the map above
    }
    dir <- file.path(OUTPUT_ROOT, folder_name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file.path(dir, filename)
}

# Builds a path like "1_Output/4_clusters/human_Mfuzz_Final_Clusters_Plot.png".
# Unlike out_path(), every dataset shares this ONE folder rather than getting
# its own -- used only for the final clustering deliverables (the final
# cluster plot, final centroid plot, and the soft clustering results CSV),
# so all three datasets' finished results sit together for easy comparison
# instead of each being buried in its own per-dataset diagnostics folder.
# The dataset name is prepended to the filename (since the folder is
# shared) so files from different datasets never collide.
final_output_path <- function(dataset, filename) {
    dir <- file.path(OUTPUT_ROOT, "4_clusters")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file.path(dir, paste0(dataset, "_", filename))
}
