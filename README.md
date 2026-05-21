# TREBL Analysis

This repository contains the analysis pipelines and notebooks for **TREBL** (Transcriptional Reporter by Endogenous Barcode Linking), a high-throughput assay that measures transcriptional activation domain (AD) activity over time using sequencing-based readouts.

## What is TREBL?

TREBL links activation domain (AD) sequence tiles to reporter barcodes (RPTR_BCs) during a library cloning step (Step 1). During an experiment, cells expressing different AD tiles are collected at multiple time points and sequenced. By counting how often each reporter barcode appears at each time point, we can infer the activation strength and kinetics (speed) of each AD tile.

## Experiments in this repository

### ChopTFs
A tiled library of chopped (short overlapping tiles) transcription factors from yeast. The goal is to measure activation activity and speed across thousands of AD tiles from many TFs simultaneously.

- Pipeline notebooks: [`notebooks/ChopTFs/pipeline/`](notebooks/ChopTFs/pipeline/)

### GCN4
A tiled library of the yeast transcription factor GCN4, split into multiple sequencing pools (Pool A, B, C). Includes both PolyT and UMI-based sequencing strategies. The analysis is analogous to ChopTFs but structured around GCN4-specific pools.

- Pipeline notebooks: [`notebooks/GCN4/pipeline_analysis/`](notebooks/GCN4/pipeline_analysis/)

## Repository Structure

```
TREBL-analysis/
├── notebooks/
│   ├── ChopTFs/
│   │   └── pipeline/          # Numbered pipeline notebooks for ChopTFs experiment
│   ├── GCN4/
│   │   └── pipeline_analysis/ # Numbered pipeline notebooks for GCN4 experiment
│   └── ...                    # Other exploratory notebooks
├── docs/                      # Sphinx documentation source
└── scripts/                   # Standalone analysis scripts
```

## Pipeline Overview

Both the ChopTFs and GCN4 pipelines follow the same general stages:

1. **Barcode mapping** (Steps 1–2): Map AD sequences and barcodes from library construction and experimental sequencing reads.
2. **Preprocessing** (Step 3): Quality filter raw FASTQ reads with fastp.
3. **TREBL experiment** (Step 4): Map experimental reads, deduplicate UMIs, and compute per-barcode read counts.
4. **Time normalization** (Steps 5–9): Normalize activity scores across time points using DBD controls, pool size, and read depth.
5. **Error propagation** (Steps 10–13): Quantify measurement uncertainty across biological and technical replicates.
6. **Speed fitting** (Steps 14–19): Fit mathematical models (exponential, logistic, Hill, first passage) to activity-vs-time data to extract speed parameters.
7. **Downstream analyses** (Steps 20+): Leak investigation, machine learning models (PARROT, NARDINI), GO analysis, cross-dataset comparisons.

## Dependencies

The pipelines use the `trebl_tools` Python package (internal lab package) along with standard scientific Python libraries: `pandas`, `numpy`, `scipy`, `seaborn`, `matplotlib`, `duckdb`, and `tqdm`.

Cluster jobs (heavy compute steps) are designed for the **Savio** HPC cluster at UC Berkeley.
