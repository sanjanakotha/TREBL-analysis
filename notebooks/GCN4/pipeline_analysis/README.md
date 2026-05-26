# GCN4 Pipeline Analysis

Analysis for the **GCN4 TREBL experiment**. The library was run across three sequencing pools (A, B, C), with Pool C run under two UMI strategies (PolyT and UMI-based). Run numbered files in order. 

---

## Library & Experiment Processing
The duckdb files referenced in this section are available in the duckdb folder on savio at /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb. They are not on Github because they are too large. 

- **`01_step1.py`** — Produces barcode-to-AD-tile lookup table used in all downstream steps.
- **`02_step2.py`** — Maps barcodes from Step 2 sequencing. Takes AD and RT FASTQs from multiple experimental runs and links them back to the Step 1 map.
- **`03a_yeast_pool_A.py`** — Maps experimental FASTQs to the Step 1 barcode map, performs UMI deduplication, and outputs per-barcode count tables for Pool A.
- **`03b_yeast_pool_A_activity_scores.ipynb`** — Computes activity scores per barcode for Pool A across time points.
- **`04a_yeast_pool_B.py`** — Same as 03a for Pool B.
- **`04b_yeast_pool_B_activity_scores.ipynb`** — Computes activity scores per barcode for Pool B across time points.
- **`05a_yeast_pool_C_polyt.py`** — Same as 03a for Pool C (PolyT capture strategy).
- **`05b_yeast_pool_C_polyt_activity_scores.ipynb`** — Computes activity scores per barcode for Pool C (PolyT).
- **`06a_yeast_pool_C_umi.py`** — Same as 03a for Pool C (UMI-based strategy). Includes logic to skip already-processed samples for efficient re-runs.
- **`06b_yeast_pool_C_umi_activity_scores.ipynb`** — Computes activity scores per barcode for Pool C (UMI-based). Uses DuckDB for scalable processing.
- **`07_activity_score_comparison.ipynb`** — Compares activity scores across all pools to assess reproducibility; merges all results into a single table for downstream use.

---

## Normalization

- **`08_DBD_time_norm.ipynb`** — Normalizes activity using DBD controls to remove systematic time-dependent drift across all pools.
- **`09_pool_size_read_depth_norm.ipynb`** — Corrects for differences in pool size (distinct barcodes) and sequencing depth across replicates, pools, and time points.

---

## QC & Error Estimation

- **`10_error_propagation.ipynb`** — Estimates variance from biological replicates (same AD tile, different barcode combinations).
- **`11_bootstrap_error_propagation.ipynb`** — Model-free bootstrap confidence intervals on activity scores. Also calculates a p-value for whether a tile is significantly active above the negative control baseline.

---

## Speed Modeling

- **`12_speed_fit_prep.ipynb`** — Filters tiles with sufficient data quality and coverage across time points; outputs `activities_to_fit.csv` used in Steps 13–16.
- **`13_speed_exponential_fit.ipynb`** — Fits an exponential model to each tile's activity-vs-time curve. **Primary speed metric used downstream.**
- **`14_speed_logistic_fit.ipynb`** — Fits a logistic (sigmoidal) model as an alternative.
- **`15_speed_hill_fit.ipynb`** — Fits a Hill function model for cooperative activation dynamics.
- **`16_speed_first_passage.ipynb`** — Estimates speed as the first time point a tile crosses an activity threshold.
- **`17_speed_fit_comparison.ipynb`** — Compares all four speed models; selects best-performing.
- **`18_leak_investigation.ipynb`** — Checks whether inactive tiles show artifactual activity and whether this affects speed estimates.

---

## Downstream Analysis

- **`19a_exponential_parrot.sh`** — Submits PARROT (recurrent neural network) training job on Savio using exponential speed data.
- **`19b_exponential_parrot.ipynb`** — Compares PARROT-predicted speed from sequence features to measured values.
- **`20_composition_speed_prediction.ipynb`** — Builds amino acid composition-based regression models for predicting AD speed and strength; compares to PARROT as an interpretable baseline.
- **`21a_nardini_on_active.py`** — Runs NARDINI (localCIDER) to compute disorder/biophysical parameters for each active AD tile in parallel on Savio.
- **`21b_submit_nardini.sh`** — Slurm submission script for the above.
- **`21d_nardini_speed_prediction.ipynb`** — Correlates NARDINI biophysical parameters (charge, hydrophobicity, IDP metrics) with measured speed to find interpretable predictors.
- **`22_activities_per_TF.ipynb`** — Summarizes activity and speed distributions across GCN4 orthologs/full-length sequences.
- **`23_merge_with_ChopTFs.ipynb`** — Merges GCN4 results with the ChopTFs dataset for cross-dataset comparisons.

---

## Helper Files

- **`gcn4_analysis_savio_job.sh`** — Main Savio job submission script
- **`plotting_helpers.py`** — Shared plotting utilities for activity score notebooks