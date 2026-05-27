# TREBL Analysis

Analysis pipelines and notebooks for TREBL.

## Getting Started

Clone this repo. The main analysis notebooks are in [`notebooks/ChopTFs/pipeline/`](notebooks/ChopTFs/pipeline/) and [`notebooks/GCN4/pipeline_analysis/`](notebooks/GCN4/pipeline_analysis/). Most of the files needed to run the analysis are included in the data and output folders. The DuckDB files are too large for GitHub, so those can be found on Savio at `/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb`; you'll need access to Savio to work with those. See the [trebl_tools docs](https://trebl-tools.readthedocs.io/en/latest/index.html) for the most up to date details on installing and using the package.

### On Savio
The `trebl_tools` conda environment is already installed in the shared project directory. Activate it and register the Jupyter kernel once before running any notebooks.

```bash
conda activate /global/scratch/projects/fc_mvslab/conda/trebl_tools

# register Jupyter kernel (run once)
python -m ipykernel install --user --name trebl_tools_shared --display-name "trebl_tools (shared)"
```

### Locally
Install `trebl_tools` from the latest release by cloning the repo, creating the conda environment from the provided YAML, and installing the package.

```bash
# clone the latest release
git clone --branch v0.1.5 --depth 1 https://github.com/staller-lab/trebl_tools.git
cd trebl_tools

# create and activate conda env
conda env create -f trebl_tools_env.yaml
conda activate trebl_tools_env

# install the package
pip install .

# register Jupyter kernel
python -m ipykernel install --user --name trebl_tools_env --display-name "trebl_tools (v0.1.5)"
```

---

## Experiments

### ChopTFs
A tiled library of short overlapping peptides covering transcription factors from yeast, designed to map activation domain activity and speed at high resolution across many TFs.
- Pipeline: [`notebooks/ChopTFs/pipeline/`](notebooks/ChopTFs/pipeline/)

### GCN4
A tiled library of the yeast transcription factor GCN4, split across three sequencing pools (A, B, C) with both PolyT and UMI-based sequencing strategies. The analysis follows the same structure as ChopTFs but is organized around the GCN4-specific pooling design.
- Pipeline: [`notebooks/GCN4/pipeline_analysis/`](notebooks/GCN4/pipeline_analysis/)

---

## Repository Structure

```
TREBL-analysis/
├── notebooks/
│   ├── ChopTFs/
│   │   └── pipeline/          # Numbered pipeline notebooks for ChopTFs
│   ├── GCN4/
│   │   └── pipeline_analysis/ # Numbered pipeline notebooks for GCN4
├── docs/                      # Sphinx documentation source
└── scripts/                   # Standalone analysis scripts
```

---

## Pipeline Overview

Both pipelines follow the same general stages:

1. **Barcode mapping** — Map AD sequences and barcodes from library construction and experimental sequencing reads.
2. **Preprocessing** — Quality filter raw FASTQs with `fastp`.
3. **TREBL experiment** — Map experimental reads, deduplicate UMIs, and compute per-barcode read counts.
4. **Time normalization** — Normalize activity scores across time points using DBD controls, empty AD constructs, pool size, and read depth.
5. **Error propagation** — Quantify measurement uncertainty across biological and technical replicates using both analytical and bootstrap approaches.
6. **Speed fitting** — Fit exponential, logistic, Hill, and first-passage models to activity-vs-time curves to extract kinetic parameters.
7. **Downstream analyses** — Leak investigation, ML-based speed prediction (PARROT, NARDINI), GO enrichment analysis, and cross-dataset comparisons.

---

## Dependencies

The pipelines use the `trebl_tools` package (see Getting Started above) along with standard scientific Python libraries: `pandas`, `numpy`, `scipy`, `seaborn`, `matplotlib`, `duckdb`, and `tqdm`.