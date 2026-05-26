# ChopTFs Pipeline

Analysis for the **ChopTFs TREBL experiment** in numbered order. 

---

## Library & Experiment Processing
The duckdb files referenced in this section are available in the duckdb folder on savio at /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb. They are not on Github because they are too large. 

- **`01_step1.py`** — Produces barcode-to-AD-tile lookup table used in all downstream steps.
- **`02_step2.py`** — *(unused here, empty. Step 2 analysis conducted at notebooks/ChopTFs/01162025 ChopTFs analysis redone with pipeline scripts.ipynb)*
- **`03a_preprocess_trebl_experiment.py`** — Runs `fastp` quality filtering on single-end reads for AD and RT reads.
- **`03b_preprocess_trebl_experiment_loss.ipynb`** — Visualizes read loss from fastp filtering for QC.
- **`04a_trebl_experiment.py`** — Core processing step: maps experimental FASTQs to the Step 1 barcode map, performs UMI deduplication, and outputs per-barcode UMI count tables per replicate and time point.
- **`04b_trebl_experiment_loss.ipynb`** — Visualizes read and barcode loss through the mapping pipeline.
- **`04c_trebl_experiment_activities.ipynb`** — Computes activity scores per barcode from UMI count tables across time points.
- **`04d_reads_per_umi.ipynb`** — QC: checks reads-per-UMI distributions.

---

## Normalization

- **`05_all_data_z-score_time_normalization.ipynb`** — Explores activity score distributions over time and applies z-score normalization to correct for global population shifts.
- **`06_DBD_time_normalization.ipynb`** — Normalizes activity using DBDs to remove systematic time-dependent drift.
- **`07a_empty_AD_time_normalization_step1.py`** — Re-runs Step 1 barcode mapping for empty AD constructs (no insert), which serve as an additional negative control.
- **`07b_empty_AD_time_normalization.ipynb`** — Incorporates empty AD counts into time normalization as a baseline signal.
- **`08_pool_size_read_depth_normalization.ipynb`** — Corrects for differences in pool size (distinct barcodes) and sequencing depth across replicates and time points.
- **`09_all_time_normalization.ipynb`** — Applies DBD normalization, empty AD correction, and read depth normalization; produces the primary activity table used in all downstream steps.

---

## QC & Error Estimation

- **`10a_downsampling_split_files.ipynb`** — Splits FASTQs into subsampled chunks at varying depths.
- **`10b_downsampling.py`** — Processes each downsampled FASTQ through the TREBL pipeline on Savio as a parallelized array job.
- **`10c_downsampling.sh`** — Slurm submission script for the above.
- **`10d_downsampling_analysis.ipynb`** — Analyzes how activity estimates change with sequencing depth; checks saturation and barcode recovery.
- **`11_error_propagation_bio_only.ipynb`** — Estimates variance from biological replicates (same AD tile, different barcode combinations).
- **`12_error_propogation_bio_and_tech.ipynb`** — Extends error propagation to include technical replicates (same barcode, multiple sequencing runs).
- **`13_error_propagation_bootstrapping.ipynb`** — Model-free bootstrap confidence intervals on activity scores.

---

## Speed Modeling

- **`14_speed_fit_prep.ipynb`** — Filters tiles with sufficient data quality and coverage across time points; outputs `activities_to_fit.csv` used in Steps 15–18.
- **`15_speed_exponential_fit.ipynb`** — Fits an exponential model to each tile's activity-vs-time curve. **Primary speed metric used downstream.**
- **`16_speed_logistic_fit.ipynb`** — Fits a logistic (sigmoidal) model as an alternative.
- **`17_speed_hill_fit.ipynb`** — Fits a Hill function model for cooperative activation dynamics.
- **`18_speed_first_passage.ipynb`** — Estimates speed as the first time point a tile crosses an activity threshold.
- **`19_speed_fit_comparison.ipynb`** — Compares all four speed models; selects best-performing.
- **`20_leak_investigation.ipynb`** — Checks whether inactive tiles show artifactual activity and whether this affects speed estimates.

---

## Downstream Analysis

- **`21a_exponential_parrot.sh`** — Submits PARROT (recurrent neural network) training job on Savio using exponential speed data.
- **`21b_exponential_parrot.ipynb`** — Compares PARROT-predicted speed from sequence features to measured values.
- **`22_comparison_to_sanborn_intervals.ipynb`** — Compares ChopTFs speed and activity to published AD annotations from Sanborn et al.
- **`23_speed_heatmap.ipynb`** — Heatmap of speed and activity across all tiles, aligned to TF sequences.
- **`24_speed_go_analysis.ipynb`** — GO enrichment analysis on TFs grouped by activation speed.
- **`25_idea_comparison.ipynb`** — Compares TREBL speed and activity to the IDEA database.

---

## Helper Files

- **`time_normalization_helpers.py`** — Shared functions for Steps 5–13
- **`speed_plotting_helpers.py`** — Shared plotting utilities for Steps 15–19
- **`chop_tf_analysis.sh`** — Savio job submission for ChopTFs analysis
- **`downsampling.sh`** — Savio job submission for downsampling