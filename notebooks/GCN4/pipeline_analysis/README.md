# GCN4 Pipeline Analysis

This directory contains the numbered pipeline notebooks and scripts for the **GCN4 TREBL experiment** — a tiled library of the yeast transcription factor GCN4, sequenced across multiple experimental pools and time points to measure activation domain (AD) activity and kinetics (speed).

The GCN4 library was run across three sequencing pools (A, B, C), with Pool C run under two different UMI strategies (PolyT and UMI-based). Run the numbered files in order. Scripts (`.py`) are intended to be run on the Savio HPC cluster; notebooks (`.ipynb`) can be run locally or on the cluster.

---

## Step-by-Step Summary

### Step 1 — Library barcode mapping
**`01_step1.py`**
Maps barcodes from the **library construction sequencing** (Step 1 FASTQ). Identifies AD sequences, AD barcodes (AD_BC), and reporter barcodes (RPTR_BC) in reads and creates a barcode-to-sequence lookup table. This map links each reporter barcode to a specific GCN4 AD tile and is used in all downstream steps.

---

### Step 2 — Experimental sequencing barcode mapping
**`02_step2.py`**
Maps barcodes from **Step 2 sequencing** (a separate library construction step for experimental reads). Takes AD and RT FASTQ files from multiple experimental runs and links them back to the Step 1 map.

---

### Step 3 — Yeast pool A TREBL experiment
**`03a_yeast_pool_A.py`** *(cluster script)*
Runs the full TREBL experiment analysis for **Yeast Pool A** (first biological replicate pool). Maps experimental FASTQ reads to the Step 1 barcode map, performs UMI deduplication, and produces per-barcode read/UMI count tables.

**`03b_yeast_pool_A_activity_scores.ipynb`**
Calculates **activity scores** per barcode combination for Pool A. Merges AD and RT counts and computes the per-tile activity metric at each time point.

---

### Step 4 — Yeast pool B TREBL experiment
**`04a_yeast_pool_B.py`** *(cluster script)*
Runs the full TREBL experiment analysis for **Yeast Pool B** (second biological replicate pool), analogous to Step 3a.

**`04b_yeast_pool_B_activity_scores.ipynb`**
Calculates activity scores per barcode combination for Pool B.

---

### Step 5 — Yeast pool C (PolyT) TREBL experiment
**`05a_yeast_pool_C_polyt.py`** *(cluster script)*
Runs the TREBL experiment analysis for **Yeast Pool C using PolyT-based sequencing** (third pool, PolyT capture strategy). Analogous to Steps 3a and 4a.

**`05b_yeast_pool_C_polyt_activity_scores.ipynb`**
Calculates activity scores per barcode combination for Pool C (PolyT).

---

### Step 6 — Yeast pool C (UMI-based) TREBL experiment
**`06a_yeast_pool_C_umi.py`** *(cluster script)*
Runs the TREBL experiment analysis for **Yeast Pool C using a UMI-based sequencing strategy** (same biological pool as Step 5, but with explicit UMI barcodes in reads). Includes logic to skip already-processed samples for efficient re-runs.

**`06b_yeast_pool_C_umi_activity_scores.ipynb`**
Calculates activity scores per barcode combination for Pool C (UMI-based). Uses a DuckDB database for scalable processing of the larger UMI dataset.

---

### Step 7 — Activity score comparison across pools
**`07_activity_score_comparison.ipynb`**
Compares activity scores across Pool A, Pool B, and Pool C (PolyT and UMI) to assess reproducibility and identify discrepancies. Merges all pool results into a single data table for downstream use.

---

### Step 8 — DBD time normalization
**`08_DBD_time_norm.ipynb`**
Normalizes activity scores using **DNA-Binding Domain (DBD) controls** — sequences in the library expected to be inactive. Removes systematic time-dependent drift across all pools, yielding normalized activity scores per tile per time point.

---

### Step 9 — Pool size and read depth normalization
**`09_pool_size_read_depth_norm.ipynb`**
Corrects for differences in pool size (number of distinct barcodes present) and sequencing depth across replicates, pools, and time points.

---

### Step 10 — Error propagation (biological replicates)
**`10_error_propagation.ipynb`**
Estimates measurement uncertainty by comparing activity scores across **biological replicates** (same AD tile, different barcode combinations). Computes replicate-level variance and standard error.

---

### Step 11 — Error propagation via bootstrapping
**`11_bootstrap_error_propagation.ipynb`**
Uses **bootstrapping** (resampling with replacement) to estimate confidence intervals on activity scores. Also calculates a p-value for whether a tile is significantly active above the negative control baseline.

---

### Step 12 — Speed fit preparation
**`12_speed_fit_prep.ipynb`**
Filters and formats the normalized activity-vs-time data for kinetic model fitting. Keeps only tiles with sufficient data quality and coverage across time points. Outputs the `activities_to_fit.csv` table used in Steps 13–16.

---

### Step 13 — Speed: exponential fit
**`13_speed_exponential_fit.ipynb`**
Fits an **exponential model** to each tile's activity-vs-time curve to extract kinetic parameters (half-time, plateau, rate constant). The primary speed metric used in downstream analyses.

---

### Step 14 — Speed: logistic fit
**`14_speed_logistic_fit.ipynb`**
Fits a **logistic (sigmoidal) model** to the activity-vs-time curves as an alternative to the exponential.

---

### Step 15 — Speed: Hill function fit
**`15_speed_hill_fit.ipynb`**
Fits a **Hill function** model to the activity-vs-time curves (cooperative activation dynamics).

---

### Step 16 — Speed: first passage time
**`16_speed_first_passage.ipynb`**
Estimates speed using a **first passage time** framework — the time point at which a tile first crosses an activity threshold.

---

### Step 17 — Speed model comparison
**`17_speed_fit_comparison.ipynb`**
Compares speed estimates from the exponential, logistic, Hill, and first-passage models. Assesses agreement across models to determine which is most appropriate.

---

### Step 18 — Leak investigation
**`18_leak_investigation.ipynb`**
Investigates **background/leak signal** — whether inactive (non-activating) tiles show artifactual activity and how this affects speed estimates for the GCN4 dataset.

---

### Step 19 — PARROT machine learning model
**`19a_exponential_parrot.sh`** *(cluster script)*
Submits a job to train/run **PARROT** (a recurrent neural network for predicting sequence-to-function relationships) on the exponential speed data on the Savio cluster. Also includes variant logs for merged ChopTFs+GCN4 runs.

**`19b_exponential_parrot.ipynb`**
Analyzes PARROT predictions — compares predicted speed from sequence features to measured speed values and evaluates model performance.

---

### Step 20 — Composition-based speed prediction
**`20_composition_speed_prediction.ipynb`**
Builds and evaluates simpler **amino acid composition-based regression models** for predicting AD speed and strength from GCN4 tile sequences. Compares these to the PARROT model as an interpretable baseline.

---

### Step 21 — NARDINI disorder-based speed prediction
**`21a_nardini_on_active.py`** *(cluster script)*
Runs the **NARDINI** tool (from localCIDER) to compute disorder/biophysical sequence parameters for each active GCN4 AD tile. Runs in parallel using all available CPUs on the Savio cluster.

**`21b_submit_nardini.sh`** *(Slurm job submission script)*
Submits the NARDINI job to the Savio cluster.

**`21d_nardini_speed_prediction.ipynb`**
Analyzes NARDINI output — correlates computed biophysical sequence parameters (e.g., charge, hydrophobicity, IDP metrics) with measured speed to find interpretable predictors.

---

### Step 22 — Activities per TF
**`22_activities_per_TF.ipynb`**
Summarizes and visualizes **activity and speed distributions per GCN4 ortholog/full-length sequence**, showing which regions of GCN4 are most active and how speed varies across the protein.

---

### Step 23 — Merge with ChopTFs data
**`23_merge_with_ChopTFs.ipynb`**
Merges GCN4 TREBL results with the **ChopTFs dataset** to enable cross-dataset comparisons of speed and activity predictions (e.g., using shared PARROT models).

---

## Helper Files

| File | Description |
|------|-------------|
| `gcn4_analysis_savio_job.sh` | Main Savio cluster job submission script |
| `plotting_helpers.py` | Shared plotting utilities (used by activity score notebooks) |
| Log files (`.log`, `.out`) | Cluster job output logs for debugging |
