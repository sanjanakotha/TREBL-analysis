# ChopTFs Pipeline

This directory contains the numbered pipeline notebooks and scripts for the **ChopTFs TREBL experiment** — a tiled transcription factor library designed to measure activation domain (AD) activity and speed across thousands of short peptide tiles from many yeast TFs.

Run the numbered files in order. Scripts (`.py`) are intended to be run on the Savio HPC cluster; notebooks (`.ipynb`) can be run locally or on the cluster.

---

## Step-by-Step Summary

### Step 1 — Library barcode mapping
**`01_step1.py`**
Maps barcodes from the **library construction sequencing** (Step 1 FASTQ). Identifies AD sequences, AD barcodes (AD_BC), and reporter barcodes (RPTR_BC) in reads and creates a barcode-to-sequence lookup table. This map links each reporter barcode to a specific AD tile and is used in all downstream steps.

---

### Step 2 — (placeholder / not used for ChopTFs)
**`02_step2.py`** *(empty — Step 2 sequencing was not needed for this experiment)*

---

### Step 3 — Preprocess experimental FASTQ reads
**`03a_preprocess_trebl_experiment.py`** *(cluster script)*
Runs `fastp` quality filtering on the raw experimental FASTQ files for both AD and RT (reporter) reads. Writes filtered reads and a summary log.

**`03b_preprocess_trebl_experiment_loss.ipynb`**
Visualizes read loss from fastp filtering — shows how many reads were removed per replicate and time point for QC.

---

### Step 4 — TREBL experiment analysis
**`04a_trebl_experiment.py`** *(cluster script)*
The core sequencing processing step. Maps experimental FASTQ reads (post-fastp) to the Step 1 barcode map, performs UMI deduplication, and produces per-barcode UMI count tables for each replicate and time point.

**`04b_trebl_experiment_loss.ipynb`**
Visualizes read and barcode loss through the TREBL experiment mapping pipeline (how many reads map, how many barcodes are recovered, etc.).

**`04c_trebl_experiment_activities.ipynb`**
Calculates **activity scores** per barcode combination from the UMI count tables. Merges AD and RT counts and computes the activity metric for each AD tile at each time point.

**`04d_reads_per_umi.ipynb`**
Quality control: inspects the distribution of reads per UMI to check for PCR over-amplification or other artifacts.

---

### Step 5 — Z-score time normalization (all data)
**`05_all_data_z-score_time_normalization.ipynb`**
Explores how the distribution of activity scores changes over time. Applies a z-score transformation across time points to normalize for global shifts in the population.

---

### Step 6 — DBD time normalization
**`06_DBD_time_normalization.ipynb`**
Normalizes activity scores using **DNA-Binding Domain (DBD) controls** — sequences in the library that are expected to be inactive. This removes systematic time-dependent drift, yielding a normalized activity score per tile per time point.

---

### Step 7 — Empty AD time normalization
**`07a_empty_AD_time_normalization_step1.py`** *(cluster script)*
Re-runs Step 1 barcode mapping specifically for **empty AD constructs** (no activation domain insert). These serve as an additional negative control for normalization.

**`07b_empty_AD_time_normalization.ipynb`**
Analyzes the empty AD read counts and incorporates them into the time normalization as a baseline signal.

---

### Step 8 — Pool size and read depth normalization
**`08_pool_size_read_depth_normalization.ipynb`**
Corrects for differences in pool size (number of distinct barcodes present) and sequencing depth across replicates and time points.

---

### Step 9 — Final combined time normalization
**`09_all_time_normalization.ipynb`**
Summarizes and visualizes the final normalized activity scores across all time points after applying DBD normalization, empty AD correction, and read depth normalization. Produces the primary activity table used downstream.

---

### Step 10 — Downsampling
**`10a_downsampling_split_files.ipynb`**
Splits FASTQ files into subsampled chunks at different depths, in preparation for downsampling analysis.

**`10b_downsampling.py`** *(cluster script)*
Processes each downsampled FASTQ file through the TREBL pipeline (barcode mapping, UMI deduplication) on the Savio cluster as a parallelized array job.

**`10c_downsampling.sh`** *(Slurm job submission script)*
Submits the downsampling jobs to the Savio cluster.

**`10d_downsampling_analysis.ipynb`**
Analyzes how activity estimates change as a function of sequencing depth — checks whether the experiment is saturating and how reliably barcodes are recovered at lower depths.

---

### Step 11 — Error propagation (biological replicates only)
**`11_error_propagation_bio_only.ipynb`**
Estimates measurement uncertainty by comparing activity scores across **biological replicates** (same AD tile, different barcode combinations). Computes replicate-level variance.

---

### Step 12 — Error propagation (biological + technical replicates)
**`12_error_propogation_bio_and_tech.ipynb`**
Extends error propagation to also account for **technical replicates** (same barcode, multiple sequencing runs). Combines both sources of variance.

---

### Step 13 — Error propagation via bootstrapping
**`13_error_propagation_bootstrapping.ipynb`**
Uses **bootstrapping** (resampling with replacement) to estimate confidence intervals on activity scores in a model-free way.

---

### Step 14 — Speed fit preparation
**`14_speed_fit_prep.ipynb`**
Filters and formats activity-vs-time data for speed model fitting. Keeps only tiles with sufficient data quality and coverage across time points. Outputs the `activities_to_fit.csv` table used in Steps 15–18.

---

### Step 15 — Speed: exponential fit
**`15_speed_exponential_fit.ipynb`**
Fits an **exponential model** to each tile's activity-vs-time curve to extract kinetic parameters (e.g., half-time, plateau). The primary speed metric used in downstream analyses.

---

### Step 16 — Speed: logistic fit
**`16_speed_logistic_fit.ipynb`**
Fits a **logistic (sigmoidal) model** to the activity-vs-time curves as an alternative to the exponential.

---

### Step 17 — Speed: Hill function fit
**`17_speed_hill_fit.ipynb`**
Fits a **Hill function** model to the activity-vs-time curves (cooperative activation dynamics).

---

### Step 18 — Speed: first passage time
**`18_speed_first_passage.ipynb`**
Estimates speed using a **first passage time** framework — the time point at which a tile first crosses an activity threshold.

---

### Step 19 — Speed model comparison
**`19_speed_fit_comparison.ipynb`**
Compares the speed estimates from the exponential, logistic, Hill, and first-passage models. Assesses agreement and selects the best-performing model.

---

### Step 20 — Leak investigation
**`20_leak_investigation.ipynb`**
Investigates **background/leak signal** — whether inactive (non-activating) tiles show artifactual activity and whether this affects speed estimates.

---

### Step 21 — PARROT machine learning model
**`21a_exponential_parrot.sh`** *(cluster script)*
Submits a job to train/run **PARROT** (a recurrent neural network for predicting properties from protein sequences) on the exponential speed data on the Savio cluster.

**`21b_exponential_parrot.ipynb`**
Analyzes PARROT predictions — compares predicted speed from sequence features to measured speed values.

---

### Step 22 — Comparison to Sanborn et al. intervals
**`22_comparison_to_sanborn_intervals.ipynb`**
Compares ChopTFs speed and activity measurements to published activation domain annotations from [Sanborn et al.](https://doi.org/10.1126/science.add5701).

---

### Step 23 — Speed heatmap
**`23_speed_heatmap.ipynb`**
Visualizes speed and activity across all ChopTFs tiles as a heatmap aligned to TF sequences.

---

### Step 24 — Speed GO analysis
**`24_speed_go_analysis.ipynb`**
Performs **Gene Ontology (GO) enrichment analysis** on TFs grouped by activation speed to identify functional categories associated with fast or slow activators.

---

### Step 25 — IDEA comparison
**`25_idea_comparison.ipynb`**
Compares TREBL speed and activity data to the **IDEA database** (intrinsically disordered enhancer-like activation domains).

---

## Helper Files

| File | Description |
|------|-------------|
| `time_normalization_helpers.py` | Shared functions for time normalization (used in Steps 5–13) |
| `speed_plotting_helpers.py` | Shared plotting utilities for speed fit visualization (Steps 15–19) |
| `chop_tf_analysis.sh` | Shell script for running ChopTFs analysis jobs on Savio |
| `downsampling.sh` | Shell script for downsampling jobs |
