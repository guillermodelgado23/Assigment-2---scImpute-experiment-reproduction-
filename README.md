# Assigment-2---scImpute-experiment-reproduction-
# scImpute — revised Task 2 experiment

This project uses the **original R implementation by Li & Li** and changes **only our simulation and experimental comparison**, not the statistical scImpute algorithm.

1. We simulate three types, 20,000 genes, 150 cells, and 810 differential genes for the full run.
2. **Adaptation for realistic count input:** apply 2–10-fold changes to *counts*, rather than multiplying the log10 means (equivalent to adding `log10(fold)` on log scale). This avoids enormous simulated counts and eliminates artificial upper clipping.
3. **Adaptation for a smaller study:** use 70% of the paper-inspired dropout probability, defined here as `0.70 * exp(-0.1 * gene_mean^2)`. This is **not** the exact simulated dropout mechanism used in the publication.
4. Run scImpute in **unlabeled mode** (the primary analysis, without giving known biological cell labels).
5. Run scImpute a second time in **labeled mode** (an expressly separate *diagnostic*: gives the known labels to test the effect of the package's internal clustering). This second mode must NOT be described as a reproduction of the authors' primary unlabeled setting.
6. Calculate PCA WCSS, normalized WCSS (WCSS/TSS), and masked-entry RMSE. The diagnostics enable comparison without pretending numeric equivalence to Figure 4a.

## Run in RStudio (Windows)

1. Save `scImpute experiment simulation - Guillermo Delgado.R` in your assignment folder.
2. Open RStudio and select **Session > Set Working Directory > Choose Directory** and select that folder, or run `setwd("C:/.../Assignment_2")` in the Console.
3. Open the R file and ensure:

```r
DATA_SOURCE <- "new"
QUICK_TEST <- TRUE
RUN_LABELED_DIAGNOSTIC <- TRUE
```

4. Click **Source** to execute the entire script, not selected lines.
5. Once the pilot finishes, change **only** `QUICK_TEST <- FALSE` and click **Source** again. The script creates a new timestamped folder in `results_scImpute_v2/`, for example `full_new_20261023_150000/`.
6. The first time you run it, installation may be required. R version and `scImpute` version are saved automatically.

**Optional diagnostic on the original data, without regenerating them:** set `DATA_SOURCE <- "previous"` and `QUICK_TEST <- FALSE`, provided `results_scImpute/full/complete_counts.csv` and `raw_counts.csv` are in your working directory. Both modes run on the old exact data. This comparison is useful to isolate the effects of incorrect internal clusters from the simulation problems.

## Files to upload for analysis

From the newly generated **full** directory:
- `RESULTS_SUMMARY.csv`
- `FIGURE_PCA.png`
- `CLUSTER_TABLE.csv` (produced if the author's package saves its `clust.rds`)
- `WITHIN_TYPE_WCSS.csv`
- `PARAMETERS_USED.csv`
- `REPORT_VALUES.txt`

No need to upload large 20,000 × 150 count matrices unless the new diagnostics reveal another issue.

## Academic integrity of the comparison

Li & Li (2018), *Nature Communications* 9:997, Figure 4a reported PCA WCSS values 94 (Complete), 2646 (Raw), and 152 (scImpute). Those exact values are **not** expected from this newly adapted count simulation. The revised report should present our measured WCSS values alongside the paper's values, but clearly explain that the simulation process, dropout prevalence, random seed, R package version, and possibly PCA preprocessing differ. If labeled diagnostic performs better than unlabeled, do not describe the diagnostic run as the original paper experiment.

Official implementation: https://github.com/Vivianstats/scImpute

Documentation for labeled mode: https://rdrr.io/github/Vivianstats/scImpute/f/vignettes/scImpute-vignette.Rmd
