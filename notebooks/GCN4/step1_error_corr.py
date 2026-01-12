import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
from concurrent.futures import ThreadPoolExecutor, as_completed
import dask.dataframe as dd
import sys
import glob
sys.path.append("../../scripts")
import initial_map
import map_refiner
import complexity
import finder
import preprocess
import complexity
import plotting

import glob
from tqdm import tqdm  # progress bar


print("Libraries and custom modules imported successfully.")

# Define barcode structures used to parse sequencing reads

print("Defining barcode configurations...")

EC_AD = finder.Barcode(
    name="AD",
    preceder="GGCTAGC",
    post="",
    length=120
)

EC_AD_BC = finder.Barcode(
    name="AD_BC",
    preceder="CGCGCC",
    post="",
    length=11
)

EC_RPTR_BC = finder.Barcode(
    name="RPTR_BC",
    preceder="CTCGAG",
    post="",
    length=14
)

print("Barcode definitions complete.")

# 1. Initial mapping of sequencing reads to barcodes/designs

print("Initializing Step 1 mapper...")

step1_mapper = initial_map.InitialMapper(
    db_path="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_err_corr.db",
    step_name="step1",  # Include time point / replicate in name if needed
    seq_file="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A10_S1.fastq.gz.assembled.fastq",
    design_file_path="/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/A10_sequencing/v2/current/a10_designfile.csv",
    bc_objects=[EC_AD, EC_AD_BC, EC_RPTR_BC],
    reverse_complement=True
)

print("Running initial mapping...")
step1_mapper.create_map()
print("Initial mapping completed.")

# 2. Refine mapping results with filtering and validation

print("Initializing map refiner...")

refiner = map_refiner.MapRefiner(
    db_path="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_err_corr.db",
    bc_objects=[EC_AD, EC_AD_BC, EC_RPTR_BC],
    column_pairs=[("RPTR_BC", "AD")],
    design_check=True,
    reads_threshold=50,
    map_order=[
        "barcode_exists",
        "quality",
        "error_corrected",
        "grouped",
        "thresholded",
        "unique_target",
        "designed"
    ],
    step_name="step1",
    should_check_exists=False,
    plot_histograms=True,
    output_figures_path="../../output/GCN4_error_corrected/figures/",
    design_file = "/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/A10_sequencing/v2/current/a10_designfile.csv"
)

print("Refining map from database...")
refiner.refine_map_from_db()
print("Map refinement completed.")