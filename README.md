# Comparative Time-Course Expression Clustering: Human, Mouse, and Control

RNA-seq time-course analysis comparing human, mouse, and human-control samples
via ortholog-matched gene sets, DESeq2 normalization, quality control, and
Mfuzz soft clustering. Two analyses are included: the full timecourse, and a
parallel, fully independent rerun restricted to the first 120 minutes (see
`4_first_2_hours_analysis/`).

See `5_docs/Methods_Report.odt` for the full write-up (methods, rationale,
and results, including both analyses). This README covers repo layout and
how to reproduce the pipeline.

## Repository layout

```
1_data/matched/                Gene-set-matched input CSVs (one per dataset)
                                -- the actual starting point for the
                                diagnostics scripts

2_scripts_includes_legacy/     Master copies of the current pipeline scripts
  filter_orthologs.R             Active -- filters raw counts to 1:1 orthologs
  pipeline_utils.R               Shared helper: output path/folder handling
  01_diagnostics.R                Master copy, full timecourse
  02_clustering.R                 Master copy, full timecourse
  01_diagnostics_120min.R         Master copy, <=120-minute subset
  02_clustering_120min.R          Master copy, <=120-minute subset
  legacy/                        Superseded earlier versions, kept for
                                  history/reference (includes
                                  match_gene_sets.R -- see note below)

3_full_time_analysis/          Full-timecourse analysis
  scripts/                       Pinned copies of the diagnostics/clustering
                                  scripts as actually run for this analysis
  diagnostic_results/
    data_analysis/                Sample- and gene-level QC, per dataset
    cluster_analysis/              Cluster-count selection evidence, per
                                    dataset
  results/                        Final clustering outputs (plot + results
                                   CSV), one folder per dataset

4_first_2_hours_analysis/      <=120-minute subset analysis - structured
                                identically to 3_full_time_analysis/, plus
                                its own README with sample-size caveats

5_docs/                        The methods report (includes a full file
                                manifest as an appendix)
```

Within `3_full_time_analysis/` and `4_first_2_hours_analysis/`, the
`diagnostic_results/` and `results/` subfolders are further split by dataset
(`1_human/`, `2_mouse/`, `3_control/`), each holding the same set of file
types.

## Pipeline, in order

1. **`2_scripts_includes_legacy/filter_orthologs.R`** -- filters each
   species' raw counts to genes with a confirmed 1:1 ortholog, using MGI's
   `HOM_ProteinCoding.rpt` reference table (external download, not included
   here -- see link below).
2. **`2_scripts_includes_legacy/legacy/match_gene_sets.R`** -- restricts all
   three ortholog-filtered datasets to their common gene set, row-aligned.
   Produces the files under `1_data/matched/`. This step has already been
   run and its output checked in; the script itself is now kept in `legacy/`
   for reference and only needs re-running if the matched CSVs need to be
   rebuilt from scratch.
3. **`01_diagnostics.R`** (full timecourse) / **`01_diagnostics_120min.R`**
   (<=120-minute subset) -- loads the matched CSVs, excludes flagged outlier
   samples, runs DESeq2 normalization + VST, generates sample- and
   gene-level QC diagnostics, filters genes (mean > 1, then top-quartile
   CV), and evaluates candidate cluster counts. Saves per-dataset state to
   `pipeline_state.rds` for step 4.
4. **`02_clustering.R`** / **`02_clustering_120min.R`** -- reads the saved
   state, previews cluster centroids across candidate cluster counts,
   prompts for a final cluster count per dataset, and writes the final
   clustering outputs.

The master copies of steps 3-4 live in `2_scripts_includes_legacy/`; the
copies actually run to produce the checked-in results are pinned under
`3_full_time_analysis/scripts/` and `4_first_2_hours_analysis/scripts/`
respectively. `pipeline_utils.R` is a shared helper (output path/folder
handling) sourced by all of the above -- keep it in the same folder when
running a script standalone.

### Requirements

R with the following packages: `readr`, `DESeq2`, `Mfuzz`, `Biobase`.

### Running

```r
# from 3_full_time_analysis/scripts/ (or 4_first_2_hours_analysis/scripts/)
Rscript 01_diagnostics.R      # or 01_diagnostics_120min.R
Rscript 02_clustering.R       # or 02_clustering_120min.R
```

The clustering script will pause per dataset and ask you to enter a cluster
count after reviewing the centroid previews it just wrote (or press Enter to
accept the shown default). Set `prompt_for_c_final <- FALSE` inside the
script to skip the pause and use the hardcoded `c_final_choices` instead.

## Data source

Raw RNA-seq count data used in this repository originate from:

> Barry C, Schmitz MT, Argus C, Bolin JM, Probasco MD, Leng N, Duffin BM,
> Steill J, Swanson S, McIntosh BE, Stewart R, Kendziorski C, Thomson JA,
> Bacher R (2019). Automated minute scale RNA-seq of pluripotent stem cell
> differentiation reveals early divergence of human and mouse gene
> expression kinetics. *PLoS Computational Biology* 15(12): e1007543.
> https://doi.org/10.1371/journal.pcbi.1007543

Released under a Creative Commons Attribution License (CC BY), which
permits reuse and redistribution with attribution. Original data deposited
at GEO accession
[GSE129014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129014);
the original authors' own analysis code is available at
https://github.com/rhondabacher/RobotNeuralDiffPaper.

## Diagnostic results (full timecourse)

`3_full_time_analysis/diagnostic_results/` holds the QC and
cluster-count-selection evidence referenced in the report, organized into
two subfolders, each further split by dataset (`1_human/`, `2_mouse/`,
`3_control/`):

**`data_analysis/`** -- sample- and gene-level QC:
- `*_VST_Boxplot_By_Sample.pdf` -- per-sample scale check
- `*_VST_Sample_Correlation_Heatmap.pdf` -- per-sample outlier evidence
  (this is what identified the excluded samples listed in the report's
  Table 1)
- `*_Gene_CV_Histogram.pdf` -- per-dataset CV distribution behind the
  percentile-based gene filtering threshold
- `*_Gene_Mean_Expression_Histogram.pdf` -- gene-level QC sanity check

**`cluster_analysis/`** -- cluster-count selection evidence, complete for
all three datasets: centroid previews (c = 2, 3, and 4) and the cluster
overlap evaluation PDF (c = 2, 3, 4, 6, 8, 10), for human, mouse, and
control.

To regenerate these, re-run the scripts in `3_full_time_analysis/scripts/`.

## First-two-hours subset analysis

`4_first_2_hours_analysis/` contains a separate, fully independent rerun of
this entire pipeline restricted to samples at or before 120 minutes, for
all three datasets. See its own README for sample-size caveats and run
instructions.

**Status: complete.** A final cluster count of c = 3 was chosen for all
three datasets (see the report's Section 10 for the full write-up,
including a notable early transient in mouse at X.m_0/X.m_4), and the
Mfuzz run has been completed for all three -- `results/` is now populated.
For human and control, c = 3 matches the overlap-evaluation preview; for
mouse, c = 3 overrides the overlap evaluation's own recommendation (c = 2)
for consistency across datasets. Gene-assignment confidence (median
membership / % of genes above 0.7) turned out to vindicate that call for
mouse (0.999 / 94.6%) and is solid for human (0.810 / 65.9%), but for
control it's *lower* than the full-timecourse c = 3 solution that was
rejected in this same repo for being too redundant (0.719 / 52.7% here vs.
0.79 / 62% there) -- see Table 5 in the report for the full numbers and
this caveat.

## Data not included in this repository

- **Original raw count files** (`HumanRawCount.ods`, `MouseRawCount.ods`,
  `HumanControlRawCount.ods`) -- superseded by the matched CSVs in
  `1_data/matched/`, which are what the pipeline actually consumes.
- **`HOM_ProteinCoding.rpt`** -- MGI's human-mouse ortholog reference table,
  used by `filter_orthologs.R`. Download from:
  https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt

## Final cluster counts (full timecourse)

| Dataset | Cluster count (c) | Notes |
|---|---|---|
| Human   | 3 | |
| Mouse   | 3 | |
| Control | 2 (reported as final) | A 3-cluster solution showed two clusters largely overlapping (median membership 0.79); collapsing to c = 2 improved median membership to 0.91. Both the c = 2 (final) and c = 3 (superseded, retained for comparison) outputs are included under `3_full_time_analysis/results/3_control/`. |

See the methods report for full rationale, excluded-sample lists, and
per-dataset statistics.

## Final cluster counts (<=120-minute subset)

| Dataset | Cluster count (c) | Median membership | % genes > 0.7 | Notes |
|---|---|---|---|---|
| Human   | 3 | 0.810 | 65.9% | Weaker than full-timecourse human (0.974) but solid, as expected for the smaller window/sample size. |
| Mouse   | 3 | 0.999 | 94.6% | Overrides the overlap evaluation's own c = 2 recommendation, for consistency across datasets; membership confidence strongly supports the choice regardless. |
| Control | 3 | 0.719 | 52.7% | **Caveat:** weaker than the full-timecourse control's c = 3 solution (0.79 / 62%), which was itself rejected there in favor of c = 2 for being too redundant. Used here for consistency across datasets -- weigh control's cluster assignments with this in mind. |

See the methods report (Section 10, Table 5) for full rationale and the
per-dataset sample/gene counts.
