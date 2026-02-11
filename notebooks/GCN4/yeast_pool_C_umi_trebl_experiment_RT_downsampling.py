#!/usr/bin/env python3
print("PYTHON STARTED", flush=True)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob
from tqdm import tqdm
import os
import argparse

# Add TREBL scripts to path
sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/")
from scripts import initial_map, map_refiner, complexity, finder, preprocess, error_correct, plotting, umi_deduplicate

print("Done with imports", flush=True)

# Argument parsing
parser = argparse.ArgumentParser(description="Process a single FASTQ file with TREBL pipeline")
parser.add_argument("file_path", help="Path to the FASTQ file to process")
args = parser.parse_args()
file_path = args.file_path

# Paths
yeast_pool_C_umi_output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/downsampling" 
os.makedirs(yeast_pool_C_umi_output_path, exist_ok=True)
print(f"Output path: {yeast_pool_C_umi_output_path}", flush=True)

# Create per-file DB path
base_name = os.path.basename(file_path)
name_only = base_name.split('.')[0]
db_dir = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_downsampling"
os.makedirs(db_dir, exist_ok=True)
db_path = os.path.join(db_dir, f"{name_only}.db")
print(f"Using database path: {db_path}", flush=True)

# Define barcode objects
EC_RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)

# File naming
umi_path = os.path.join(yeast_pool_C_umi_output_path, f"trebl_experiment_yeast_pool_C_umi_{name_only}")
print(f"UMI output path: {umi_path}", flush=True)

# Initial mapping
umi_mapper = initial_map.InitialMapper(
    db_path=db_path,
    step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
    seq_file=file_path,
    bc_objects=[EC_RPTR_BC],
    umi_object=RT_UMI,
    reverse_complement=True,
    design_file_path=None
)
umi_mapper.create_map()

# Refinement

figures_path = os.path.join(yeast_pool_C_umi_output_path, "figures", name_only)
os.makedirs(figures_path, exist_ok=True)

refiner = map_refiner.MapRefiner(
    db_path=db_path,
    bc_objects=[EC_RPTR_BC],
    column_pairs=[],
    map_order=['quality'],
    step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
    descriptor="",
    output_figures_path=figures_path,
    reads_threshold=0,
    umi_object=RT_UMI
)
refiner.refine_map_from_db(should_check_exists=True)
refiner.plot_loss()

# Deduplication
deduplicator = umi_deduplicate.UMIDeduplicator(
    db_path=db_path,
    bc_objects=[EC_RPTR_BC],
    step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
    descriptor="",
    step1_map_name=None,
    fastq_path=file_path,
    output_path=umi_path, 
    refined_map_suffix='quality'
)
deduplicator.run_simple_deduplication()
deduplicator.save_simple_deduplication()
