# TREBL Analysis

[![Documentation Status](https://readthedocs.org/projects/trebl-analysis/badge/?version=latest)](https://trebl-analysis.readthedocs.io/en/latest/?badge=latest)

Python scripts and analysis tools for TREBL (Transcription factor REporter Barcode Library) sequencing data.

## Overview

TREBL is a high-throughput method for mapping transcription factor activation domains to their corresponding DNA barcodes through sequencing-based analysis. This repository provides a complete pipeline for processing TREBL sequencing data, including:

- **Barcode extraction and mapping**: Identify and extract barcodes from sequencing reads using customizable patterns
- **Quality control and filtering**: Apply quality thresholds, remove barcode collisions, and filter low-confidence reads
- **Error correction**: Optional UMI-based error correction to improve data quality
- **Data analysis and visualization**: Generate reports, plots, and statistics for TREBL experiments

## Features

- **Flexible barcode definitions**: Support for multiple barcode types with custom flanking sequences
- **Efficient data processing**: Uses DuckDB for fast, memory-efficient processing of large sequencing files
- **Multi-step refinement**: Configurable filtering pipeline with multiple quality control steps
- **Comprehensive documentation**: Detailed protocols and API reference available at [Read the Docs](https://trebl-analysis.readthedocs.io/)
- **Visualization tools**: Built-in plotting functions to track data loss and quality metrics

## Installation

### Prerequisites

- Python 3.10 or higher
- pip package manager

### Install Dependencies

Clone this repository and install the required packages:

```bash
git clone https://github.com/sanjanakotha/TREBL-analysis.git
cd TREBL-analysis
pip install -r docs/requirements.txt
```

### Required Packages

- pandas
- dask
- duckdb
- biopython
- matplotlib
- seaborn
- pyarrow
- tqdm
- numpy

## Quick Start

Here's a minimal example to get started with TREBL analysis using the recommended pipeline wrapper:

### Using TreblPipeline (Recommended)

```python
from scripts.pipelines import TreblPipeline
from scripts import finder

# 1. Define your barcodes
bc_objects = [
    finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
    finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
    finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
]

# 2. Initialize the pipeline
pipeline = TreblPipeline(
    db_path="./my_analysis.db",
    design_file_path="./design.csv",
    output_path="./output",  # Results saved here
)

# 3. Run complete Step 1 analysis
pipeline.run_step1(
    seq_file="./data/sequences.fastq.gz",
    bc_objects=bc_objects,
    reverse_complement=True,
    reads_threshold=50,
)

print("Analysis complete! Check ./output for results.")
```

That's it! The pipeline automatically handles:
- Barcode extraction
- Quality filtering
- Collision removal
- Result visualization
- Data export

### Alternative: Manual Step-by-Step Approach

For more control over each step:

```python
import sys
sys.path.append("./scripts")
from scripts import finder

# Define barcode structure
EC_AD = finder.Barcode(
    name="AD",
    preceder="GGCTAGC",
    post="",
    length=120,
)

EC_AD_BC = finder.Barcode(
    name="AD_BC",
    preceder="CGCGCC",
    post="",
    length=11,
)

bc_objects = [EC_AD, EC_AD_BC]
```

### 2. Run Initial Mapping

```python
from scripts import initial_map

# Create initial mapper
mapper = initial_map.InitialMapper(
    db_path="./my_analysis.db",
    step_name="step1",
    seq_file="path/to/your/sequences.fastq.gz",
    design_file_path="path/to/design_file.csv",
    bc_objects=bc_objects,
    reverse_complement=True,
)

# Run mapping (creates DuckDB database)
mapper.run()
```

### 3. Refine the Map

```python
from scripts import map_refiner

# Set up refinement parameters
refiner = map_refiner.MapRefiner(
    db_path="./my_analysis.db",
    bc_objects=bc_objects,
    step_name="step1",
    reads_threshold=50,
    column_pairs=[("AD_BC", "AD")],
    map_order=[
        "grouped",
        "thresholded",
        "barcode_exists",
        "quality",
        "unique_target",
        "designed"
    ],
)

# Run refinement
refiner.refine_map_from_db()

# Visualize results
refiner.plot_loss()
```

## Documentation

Full documentation is available at [https://trebl-analysis.readthedocs.io/](https://trebl-analysis.readthedocs.io/)

Key documentation sections:
- [Step 1 Protocol](https://trebl-analysis.readthedocs.io/en/latest/step1_protocol.html) - Detailed walkthrough for initial mapping and refinement
- [API Reference](https://trebl-analysis.readthedocs.io/en/latest/scripts.html) - Complete module and function documentation

## Project Structure

```
TREBL-analysis/
├── scripts/              # Core analysis modules
│   ├── finder.py        # Barcode definition and extraction
│   ├── initial_map.py   # Initial sequence mapping
│   ├── map_refiner.py   # Map refinement and filtering
│   ├── pipelines.py     # End-to-end pipeline workflows
│   ├── error_correct.py # Error correction utilities
│   ├── umi_deduplicate.py # UMI-based deduplication
│   ├── plotting.py      # Visualization functions
│   └── ...
├── notebooks/           # Jupyter notebooks with examples
├── docs/               # Sphinx documentation source
├── savio_jobs/         # HPC job scripts
└── README.md           # This file
```

## Usage Examples

### Using the Complete Pipeline

For more complex analyses, use the `TreblPipeline` class:

```python
from scripts.pipelines import TreblPipeline

pipeline = TreblPipeline(
    db_path="./analysis.db",
    design_file_path="design.csv",
    error_correction=True,
    output_path="./output",
)

# Run Step 1 analysis
pipeline.run_step1(
    seq_file="step1_data.fastq.gz",
    bc_objects=bc_objects,
    reverse_complement=True,
    reads_threshold=50,
)
```

### Analyzing Multiple Samples

The pipeline supports batch processing of multiple sequencing files:

```python
# Process multiple files in one run
mapper_specs = [
    {
        "seq_file": "sample1.fastq.gz",
        "step_name": "step1_sample1",
        "bc_objects": bc_objects,
        "reverse_complement": True,
        "design_file_path": "design.csv",
    },
    {
        "seq_file": "sample2.fastq.gz",
        "step_name": "step1_sample2",
        "bc_objects": bc_objects,
        "reverse_complement": True,
        "design_file_path": "design.csv",
    },
]

pipeline._run_initial_mappers(mapper_specs)
```

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

### Development Setup

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Run tests (if available)
5. Commit your changes (`git commit -am 'Add new feature'`)
6. Push to the branch (`git push origin feature/my-feature`)
7. Create a Pull Request

### Code Style

- Follow PEP 8 guidelines for Python code
- Add docstrings to all functions and classes
- Include type hints where appropriate
- Update documentation for new features

## Citation

If you use this software in your research, please cite:

```
TREBL Analysis Tools
Author: Sanjana Kotha
Year: 2025
URL: https://github.com/sanjanakotha/TREBL-analysis
```

## License

This project is available for academic and research use. Please contact the author for commercial use inquiries.

## Support

For questions, issues, or feature requests:
- Open an issue on [GitHub](https://github.com/sanjanakotha/TREBL-analysis/issues)
- Check the [documentation](https://trebl-analysis.readthedocs.io/)
- Review example notebooks in the `notebooks/` directory

## Acknowledgments

- Developed in the MVS Lab
- Built with DuckDB, pandas, and BioPython
- Documentation powered by Sphinx and Read the Docs
