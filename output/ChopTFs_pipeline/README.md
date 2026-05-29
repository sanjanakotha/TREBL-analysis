# ChopTFs Pipeline Outputs

Main output directory: `/output/ChopTFs_pipeline`

---

## TL;DR — Final Data

The primary output to use is:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_and_no_CIs.csv

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

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_active_pval.csv

For final speed and strength fits, see:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/exponential_fit.csv

---

## Table of Contents

1. [Step 1 Map](#step-1-map)
2. [Preseq Results](#preseq-results)
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

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/step1.csv

---

## Preseq Results

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/AD_transcript_size_estimate.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/RT_transcript_size_estimate.csv

### Columns

| Column              | Description                              |
|---------------------|------------------------------------------|
| `pop_size_estimate` | Estimated total transcript pool size     |
| `lower_ci`          | Lower confidence interval bound          |
| `upper_ci`          | Upper confidence interval bound          |
| `rep`               | Technical replicate                      |
| `time`              | Time in minutes                          |

---

## TREBL Experiment Raw Results

Base path:

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/trebl_experiment_/

Contains one subfolder per file. Each subfolder contains:

| File                           | Description                                                                                                       |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------|
| `preseq_input.txt`             | Read count per unique barcode (input to preseq)                                                                   |
| `preseq_yield.txt`             | Population size estimate and confidence intervals (output from preseq)                                            |
| `*_reads_per_umi.tsv`          | Read count per unique UMI using simple UMI deduplication                                                          |
| `*_simple_umi_counts.tsv`      | Simple UMI counts per barcode                                                                                     |
| `*_directional_umi_counts.tsv` | Directionally deduplicated UMI counts per barcode                                                                 |
| `*_loss_summary.csv` / `.png`  | Loss table. See [TREBL-tools docs](https://trebl-tools.readthedocs.io/en/latest/user_guide/step1.html#interpreting-outputs) for explanation of intermediates |

---

## Aggregated Results

Aggregated versions of `*_simple_umi_counts.tsv` and `*_directional_umi_counts.tsv` across all files
**(before merging with Step 1):**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/RT_trebl_experiment_results.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/AD_trebl_experiment_results.csv

> **Do not use:**
> `ChopTFs_pipeline/RPTR_BC_trebl_experiment_results.csv` — unused intermediate.

---

## Activities Per Barcode

**After merging with Step 1:**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/trebl_experiment_activities_per_barcode.csv

### Columns

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

### Paths

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/full_data_time_normalized_activities.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/empty_AD_time_normalized_activities.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/pool_size_read_depth_norm_activities.csv

### Columns

| Column               | Description                                                         |
|----------------------|---------------------------------------------------------------------|
| `mu`                 | Mean of the gaussian fit to the inactive distribution               |
| `sigma`              | Standard deviation of the gaussian fit to the inactive distribution |
| `Z-scored_activity`  | Activity Z-scored relative to the inactive distribution             |
| `shifted_activity`   | Activity shifted to the mean of the inactive distribution           |

---

## Time-Normalized Summaries

### Paths

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_fit_summaries.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/empty_AD_time_normalized_fit_summaries.csv

### Columns

| Column           | Description                                                         |
|------------------|---------------------------------------------------------------------|
| `mu`             | Mean of the gaussian fit to the inactive distribution               |
| `sigma`          | Standard deviation of the gaussian fit to the inactive distribution |
| `z-score_active` | Proportion of tiles classified as active (Z-score method)           |
| `shifted_active` | Proportion of tiles classified as active (shifted method)           |

---

## Error Propagation and Bootstrapping

### Early Attempts (unsuccessful — do not use)

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_error_propagated.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bio_error.csv

> Attempted to combine biological and technical SEs using MixedLM. Did not work when
> biological error was smaller than technical error.

### Bootstrapped Results (use these)

**Summary (mean + CI per ADseq per time):**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error.csv

**All bootstrap statistics (parquet):**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_all_bootstraps.parquet

**With classification and p-values:**

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_active_pval.csv

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

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/time_normalization/DBD_time_normalized_activities_bootstrapped_error_and_no_CIs.csv

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

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/

### Input and Aggregated Metrics

| File                  | Description                                                                          |
|-----------------------|--------------------------------------------------------------------------------------|
| `activities_to_fit.csv` | Subset of activities to which kinetic curves are fit                               |
| `aggreg_per_TF.csv`   | Speed metrics aggregated per TF. `Coverage_fraction` is the proportion of the TF covered by at least 1 tile |
| `fast_genes.txt`      | Gene IDs used as the "fast" set for GO analysis                                      |
| `background.txt`      | Background gene IDs used for GO analysis                                             |

### Baseline-Normalized Model Fits

Curves all normalized to zero at t = 0.

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/exponential_fit_baseline_norm.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/hill_fit_baseline_norm.csv
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/logistic_fit_baseline_norm.csv

### First Passage Time

| File                         | Description                                                                               |
|------------------------------|-------------------------------------------------------------------------------------------|
| `first_passage_DBD.csv`      | Uses the DBD distribution to determine when tiles first become active                     |
| `first_passage_pval_class.csv` | Uses bootstrapped differences from t=0 to determine when tiles first become active      |

### Model Fit Comparison

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/merged_speed_summary.csv

Comparison of the baseline-normalized and first passage time approaches above.

### Final Speed and Strength Fits (use this)

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/exponential_fit.csv

Fits kinetic curves of the form `f(t) = A * (1 - e^(-kt)) + C`.

| Column      | Description                                                               |
|-------------|---------------------------------------------------------------------------|
| `A`         | Amplitude parameter of the exponential curve                             |
| `k`         | Rate constant — used as the **speed** metric                             |
| `C`         | Offset parameter of the exponential curve                                |
| `A+C`       | Amplitude + offset — used as the **strength** metric                     |
| `R_squared` | Model performance relative to the data                                   |
| `flag`      | Whether the data passed the quality check                                |

### PARROT Models

    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/ChopTFs_pipeline/speed/PARROT/

Models are named by their input metric:

| Model name                    | Speed/strength definition used as input |
|-------------------------------|-----------------------------------------|
| `exponential_speed_input`     | Uses `A * k` as speed                  |
| `exponential_strength_input`  | Uses `A + C` as strength               |
| `exponential_speed_k`         | Uses `k` as speed                      |