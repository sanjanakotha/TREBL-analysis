# GCN4 Pipeline Outputs

Main output directory: `/output/GCN4_pipeline`

---

## TL;DR — Final Data

The primary output to use is:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_and_no_CIs.csv

This file contains time-normalized, bootstrapped activity measurements for each AD sequence
at each time point. Key columns:

| Column    | Description                                          |
|-----------|------------------------------------------------------|
| `ADseq`   | AD sequence identifier                               |
| `time`    | Time point                                           |
| `mean`    | Mean activity across all bootstraps                  |
| `ci_low`  | Lower bound of the confidence interval               |
| `ci_hi`   | Upper bound of the confidence interval               |

> **Note:** Rows with no `ci_low` / `ci_hi` lacked sufficient data to bootstrap and are
> lower confidence.

For classification of tiles as activators, repressors, or not significant, see:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_active_pval.csv

For final speed and strength fits, see:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/exponential_fit.csv

---

## Table of Contents

1. [Step 1 Map](#step-1-map)
2. [Step 2](#step-2)
3. [TREBL Experiment Raw Results](#trebl-experiment-raw-results)
4. [Aggregated Results](#aggregated-results)
5. [Activities Per Barcode](#activities-per-barcode)
6. [Time-Normalized Results](#time-normalized-results)
7. [Time-Normalized Summaries](#time-normalized-summaries)
8. [Error Propagation and Bootstrapping](#error-propagation-and-bootstrapping)
9. [Final Data](#final-data)
10. [Speed Fits](#speed-fits)

---

## Step 1 Map

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/step1.csv

---

## Step 2

| File                      | Description                                          |
|---------------------------|------------------------------------------------------|
| `step2_AD.csv`            | AD barcode counts                                    |
| `step2_RT.csv`            | RT barcode counts                                    |
| `step2_AD_RT_overlap.csv` | Overlap of AD and RT barcodes using the Step 1 map  |

Base path for all three:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/

---

## TREBL Experiment Raw Results

One folder per pool experiment:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_A_
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_B_
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_C_PolyT_
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_C_umi_

Each folder contains one subfolder per file.

### PolyT Experiment Subfolders

Contain visuals of read distributions and the loss table.

### UMI Experiment Subfolders

| File                           | Description                                                                                                       |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------|
| `*_reads_per_umi.tsv`          | Read count per unique UMI using simple UMI deduplication                                                          |
| `*_simple_umi_counts.tsv`      | Simple UMI counts per barcode                                                                                     |
| `*_directional_umi_counts.tsv` | Directionally deduplicated UMI counts per barcode                                                                 |
| `*_loss_summary.csv` / `.png`  | Loss table. See [TREBL-tools docs](https://trebl-tools.readthedocs.io/en/latest/user_guide/step1.html#interpreting-outputs) for explanation of intermediates |

---

## Aggregated Results

Aggregated versions of `*_simple_umi_counts.tsv` and `*_directional_umi_counts.tsv` across all
files **(before merging with Step 1).**

> **Note:** Aggregation did not work for Pool C UMI because files were too large.
> Those intermediates must be aggregated manually from the per-file subfolders.

### AD

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/AD_trebl_experiment_pool_A_results.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/AD_trebl_experiment_pool_B_results.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/AD_trebl_experiment_pool_C_PolyT_results.csv

### RT

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/RT_trebl_experiment_pool_A_results.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/RT_trebl_experiment_pool_B_results.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/RT_trebl_experiment_pool_C_PolyT_results.csv

> **Do not use:**
> `GCN4_pipeline/RPTR_BC_trebl_experiment_pool_C_umi_results.csv` — unused intermediate.

---

## Activities Per Barcode

**After merging with Step 1:**

### PolyT Pools (activity = RT / AD reads)

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_A_results_per_barcode.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_B_results_per_barcode.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_C_PolyT_results_per_barcode.csv

### UMI Pool (Pool C)

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/trebl_experiment_pool_C_umi_results_per_barcode.csv

### Columns (UMI pool)

| Column                 | Description                                          |
|------------------------|------------------------------------------------------|
| `activity_simple`      | `log10(count_simple_RT / count_simple_AD)`           |
| `activity_directional` | `log10(count_directional_RT / count_directional_AD)` |

> **Do not use:**
> `ChopTFs_pipeline/trebl_experiment_activities_per_AD.csv` — old aggregation method.
> Use the time-normalized and error-corrected version instead (see below).

---

## Time-Normalized Results

Activities Z-scored and shifted to the inactive distribution.

### Path

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities.csv

### Columns

| Column               | Description                                                         |
|----------------------|---------------------------------------------------------------------|
| `mu`                 | Mean of the gaussian fit to the inactive distribution               |
| `sigma`              | Standard deviation of the gaussian fit to the inactive distribution |
| `Z-scored_activity`  | Activity Z-scored relative to the inactive distribution             |
| `shifted_activity`   | Activity shifted to the mean of the inactive distribution           |

---

## Time-Normalized Summaries

### Path

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_fit_summaries.csv

### Columns

| Column           | Description                                                         |
|------------------|---------------------------------------------------------------------|
| `mu`             | Mean of the gaussian fit to the inactive distribution               |
| `sigma`          | Standard deviation of the gaussian fit to the inactive distribution |
| `z-score_active` | Proportion of tiles classified as active (Z-score method)           |
| `shifted_active` | Proportion of tiles classified as active (shifted method)           |

---

## Error Propagation and Bootstrapping

### Early Attempt (unsuccessful — do not use)

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities_error_propagated.csv

> Attempted to combine biological and technical SEs using MixedLM. Did not work when
> biological error was smaller than technical error.

### Bootstrapped Results (use these)

**Summary (mean + CI per ADseq per time):**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error.csv

**All bootstrap statistics (parquet):**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_all_bootstraps.parquet

**With classification and p-values:**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_active_pval.csv

### Columns (classification file)

| Column      | Description                                                                                              |
|-------------|----------------------------------------------------------------------------------------------------------|
| `class`     | `"activator"`, `"repressor"`, or `"ns"` (not significant), based on bootstrapped t=0 vs t=X differences |
| `Mean_diff` | Mean difference between t=0 and t=X                                                                     |
| `pval_fdr`  | FDR-corrected p-value                                                                                    |

> **Note:** Null values in `Mean_diff` and `pval_fdr` indicate either that the row is t=0,
> or that no data existed at t=0 and a comparison could not be made.

---

## Final Data

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_and_no_CIs.csv

Includes data even for ADseqs without sufficient barcodes to bootstrap (lower confidence).

### Columns

| Column    | Description                                                              |
|-----------|--------------------------------------------------------------------------|
| `ADseq`   | AD sequence identifier                                                   |
| `time`    | Time point                                                               |
| `mean`    | Mean activity across all bootstraps                                      |
| `ci_low`  | Lower bound of the confidence interval                                   |
| `ci_hi`   | Upper bound of the confidence interval                                   |

> **Note:** Rows with no `ci_low` / `ci_hi` did not have sufficient data points to bootstrap.
> These rows are lower confidence.

---

## Speed Fits

Base directory:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/

### Input

| File                    | Description                                       |
|-------------------------|---------------------------------------------------|
| `activities_to_fit.csv` | Subset of activities to which kinetic curves are fit |

### Baseline-Normalized Model Fits

Curves all normalized to zero at t = 0.

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/exponential_fit_baseline_norm.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/hill_fit_baseline_norm.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/logistic_fit_baseline_norm.csv

### First Passage Time

| File                          | Description                                                                          |
|-------------------------------|--------------------------------------------------------------------------------------|
| `first_passage_pval_class.csv` | Uses bootstrapped differences from t=0 to determine when tiles first become active  |

### Final Speed and Strength Fits (use this)

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/exponential_fit.csv

Fits kinetic curves of the form `f(t) = A * (1 - e^(-kt)) + C`.

| Column      | Description                                                |
|-------------|------------------------------------------------------------|
| `A`         | Amplitude parameter of the exponential curve              |
| `k`         | Rate constant — used as the **speed** metric              |
| `C`         | Offset parameter of the exponential curve                 |
| `A+C`       | Amplitude + offset — used as the **strength** metric      |
| `R_squared` | Model performance relative to the data                    |
| `flag`      | Whether the data passed the quality check                 |

### PARROT Models

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/

Models are named by their input metric and sequence set used:

| Model name                           | Description                                                               |
|--------------------------------------|---------------------------------------------------------------------------|
| `exponential_speed_input`            | Uses `A * k` as speed                                                    |
| `exponential_strength_input`         | Uses `A + C` as strength                                                 |
| `exponential_speed_k`                | Uses `k` as speed; high-confidence GCN4 sequences only                   |
| `exponential_strength_input_more_seqs` | Uses `A + C` as strength; both high- and low-confidence GCN4 sequences |
| `exponential_speed_input_k_more_seqs` | Uses `k` as speed; both high- and low-confidence GCN4 sequences        |
