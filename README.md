# TREBL Analysis

This repository contains my analysis pipelines and notebooks for TREBL.

## Experiments in this repository

### ChopTFs
A tiled library of short overlapping tiles of transcription factors from yeast.
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
2. **Preprocessing**: Quality filter raw FASTQ reads with fastp.
3. **TREBL experiment**: Map experimental reads, deduplicate UMIs, and compute per-barcode read counts.
4. **Time normalization**: Normalize activity scores across time points using DBD controls, pool size, and read depth.
5. **Error propagation**: Quantify measurement uncertainty across biological and technical replicates.
6. **Speed fitting**: Fit mathematical models (exponential, logistic, Hill, first passage) to activity-vs-time data to extract speed parameters.
7. **Downstream analyses**: Leak investigation, machine learning models (PARROT, NARDINI), GO analysis, cross-dataset comparisons.

## Dependencies

The pipelines use the `trebl_tools` Python package (internal lab package) along with standard scientific Python libraries: `pandas`, `numpy`, `scipy`, `seaborn`, `matplotlib`, `duckdb`, and `tqdm`.
