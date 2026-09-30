# First-Two-Hours Subset Analysis

A separate, fully independent analysis restricted to samples at or before
**120 minutes** (the first two hours) post-differentiation, for all three
datasets. This mirrors the main pipeline exactly (same outlier-exclusion
criteria, same normalization, QC, gene filtering, and clustering approach)
but re-runs everything from scratch on the smaller time window, rather than
slicing the full-timecourse results -- size factors, gene filtering
thresholds, and cluster assignments are all recomputed on just this subset.

## Why a separate run, not a slice

Normalization (DESeq2 size factors) and the CV-percentile gene filter are
both computed *relative to the sample set provided*. Slicing an
already-normalized, already-filtered result to <=120min would carry over
normalization and gene-selection decisions made using the full time
course, which may not be appropriate once the later timepoints are
excluded. Re-running from raw counts ensures this subset is normalized and
filtered fairly on its own terms.

## Sample sizes (read this before interpreting results)

At <=120 minutes:
- **Human**: 13 timepoints (10-minute intervals: 0, 10, 20, ..., 120)
- **Mouse**: 31 timepoints (4-minute intervals: 0, 4, 8, ..., 120)
- **Control**: 13 timepoints (10-minute intervals, same as human)

These are substantially smaller sample sizes than the full run (57/143/60
after exclusions). The CV-percentile gene filter and cluster stability 
are both more sensitive to noise with this few samples, especially for human 
and control at only 13 points each.

None of the outlier samples excluded in the main pipeline
(`NBplusNog_450/530/290/370`, `X.m_152/468`, `X.E8_410`) fall within this
window, so the same exclusion list is applied here for consistency but has
no actual effect on which samples are retained.

## Sample QC findings

No new outlier samples were identified in any of the three datasets within
this window. One finding worth flagging: mouse's two earliest samples
(`X.m_0`, `X.m_4`) show markedly reduced pairwise correlation to the rest of
the mouse dataset (mean correlation 0.895-0.904, vs. 0.95-0.98 for every
other timepoint in the window), and correlate more weakly with *each other*
(r = 0.844) than either does with any other sample. Their VST summary
statistics (mean/median/sd) are entirely normal, ruling out a global
normalization issue, and human shows no equivalent effect at its own
earliest timepoints (correlations stay tight, 0.987-0.994, throughout). This
looks like a genuine, rapid, transient response in the first four minutes of
the mouse timecourse -- consistent with the early divergence in human/mouse
kinetics described in the source paper -- rather than a technical artifact,
and `X.m_0`/`X.m_4` have been kept in rather than excluded.

## Cluster overlap evaluation

Candidate cluster counts of c = 2, 3, 4, 6, 8, and 10 were evaluated for
each dataset using the same PCA-projected centroid overlap method as the
full-timecourse analysis. Results differ markedly by dataset:

| Dataset | Well-separated at c = | First overlap detected |
|---|---|---|
| Human   | none -- overlap present at every c tested | c = 2 (smallest tested) |
| Mouse   | 2, 3, 4, 6 | c = 8 |
| Control | none -- overlap present at every c tested | c = 2 (smallest tested) |

Mouse holds up well through a wider range of cluster counts than the
full-timecourse mouse solution did. In both its c = 2 and c = 3 solutions,
the dominant feature of every cluster centroid is the large, opposite-signed
excursion at X.m_0/X.m_4 described above, converging to near-baseline
fluctuation for the rest of the window; at c = 3, two of the three clusters
differ mainly in the size of this early excursion rather than its shape and
converge to near-identical trajectories afterward, which suggests c = 2 is
the more parsimonious choice pending final confirmation.

Human and control, by contrast, show overlap already at the minimum
candidate count (c = 2), and at every larger c tested -- indicating that
within this narrower 120-minute window, genes in these two datasets don't
yet separate into distinguishable temporal groups strongly enough to
support confident multi-cluster assignment. This is itself a notable
cross-species finding, consistent with the source paper's account of
faster early transcriptional kinetics in mouse.

## Final clustering results

A final cluster count of **c = 3 was used for all three datasets**, for
consistency, and the Mfuzz run is complete:

| Dataset | c | Median membership | % genes > 0.7 | Notes |
|---|---|---|---|---|
| Human   | 3 | 0.810 | 65.9% | Solid; weaker than full-timecourse human (0.974) as expected for the smaller window/sample size. |
| Mouse   | 3 | 0.999 | 94.6% | Overrides the overlap evaluation's own c = 2 recommendation (above), but membership confidence strongly supports c = 3 regardless -- the visual similarity between Clusters 1 and 3 in the centroid preview doesn't translate into ambiguous gene assignments. |
| Control | 3 | 0.719 | 52.7% | **Caveat:** lower than the full-timecourse control's c = 3 solution (0.79 / 62%) -- the solution that was rejected there in favor of c = 2 for being too redundant. Used here anyway for cross-dataset consistency; treat control's cluster assignments accordingly. |

Genes retained (post mean>1 and CV filters): human 1,220; mouse 1,107;
control 1,131. All samples were retained (13/13 human, 31/31 mouse, 13/13
control) -- see "Sample QC findings" above.

## Running

```
cd scripts
Rscript 01_diagnostics_2hr.R
Rscript 02_clustering_2hr.R
```

Requires `../../1_data/matched/*.csv` to exist (i.e. run from inside a full
clone of the parent repo, not this folder standalone).

## Layout

```
4_first_2_hours_analysis/
├── scripts/
│   ├── 01_diagnostics_2hr.R
│   └── 02_clustering_2hr.R
├── diagnostic_results/
│   ├── data_analysis/        Populated -- CV histogram, mean-expression
│   │   ├── 1_human/           histogram, VST boxplot, correlation
│   │   ├── 2_mouse/           heatmap, per dataset (see main repo's
│   │   └── 3_control/         diagnostic_results/ for what each contains)
│   └── cluster_analysis/     Populated -- centroid previews (c = 2, 3, 4)
│       ├── 1_human/           and cluster-overlap evaluation PDF, per
│       ├── 2_mouse/           dataset
│       └── 3_control/
└── results/
    ├── 1_human/               Populated -- final cluster plot (c = 3) and
    ├── 2_mouse/                 Mfuzz soft-clustering results CSV, per
    └── 3_control/               dataset (see "Final clustering results")
```

