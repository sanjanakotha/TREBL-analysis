#!/usr/bin/env python3
print("PYTHON STARTED", flush=True)

import sys
import os
import glob
import argparse  # Necessary for argument parsing

sys.path.insert(0, "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/lib/python3.11/site-packages")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import duckdb
print(duckdb.__file__)
print(duckdb.__version__)

from tqdm import tqdm

sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/")
from scripts import initial_map, map_refiner, complexity, finder, preprocess, error_correct, plotting, umi_deduplicate, pipelines

print("Done with imports", flush=True)

# Argument parsing
parser = argparse.ArgumentParser(description="Process a single FASTQ file with TREBL pipeline")
parser.add_argument("file_path", help="Path to the FASTQ file to process")
parser.add_argument("--mode", choices=["AD", "RT"], required=True, help="Mode of operation: AD or RT")
args = parser.parse_args()
file_path = args.file_path
mode = args.mode

# Paths
yeast_pool_C_umi_output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/NKX2-2/downsampling" 
os.makedirs(yeast_pool_C_umi_output_path, exist_ok=True)
print(f"Output path: {yeast_pool_C_umi_output_path}", flush=True)


# Base filename
base_name = os.path.basename(file_path)

# Split off everything after .assembled
name_before_assembled, *rest = base_name.split(".assembled")

# Combine with the chunks/part info that comes after .assembled
chunks_part = rest[0] if rest else ""

# Clean name: remove leading/trailing dots or underscores
name_only = f"{name_before_assembled}{chunks_part}".strip("._")

db_dir = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/NKX2-2_downsampling"
os.makedirs(db_dir, exist_ok=True)
db_path = os.path.join(db_dir, f"{name_only}.db")
print(f"Using database path: {db_path}", flush=True)

# Define barcodes
trebl_exp_ADBC2 = finder.Barcode(name="ADBC2", preceder="TATGCTAT", post="GGCCGGCCG", length=6)
trebl_exp_HawkBCs = finder.Barcode(name="HawkBCs", preceder="TAGC", post="CTCGAGA", length=9)
trebl_exp_RTBC = finder.Barcode(name="RTBC", preceder="GCCCC", post="GCGG", length=16)

# UMI barcode (we just want last 12 positions so only specifying length here)
UMI = finder.Barcode(name="UMI", preceder="", post="", length=12)

# Set bc_objects based on mode
if mode == "RT":
    bc_objects = [trebl_exp_RTBC]
else:  # mode == "AD"
    bc_objects = [trebl_exp_ADBC2, trebl_exp_HawkBCs]

print(f"Mode: {mode}, Barcodes: {[bc.name for bc in bc_objects]}", flush=True)

# File naming
umi_path = os.path.join(yeast_pool_C_umi_output_path, f"trebl_experiment_{name_only}")
print(f"UMI output path: {umi_path}", flush=True)

# Initial mapping
umi_mapper = initial_map.InitialMapper(
    db_path=db_path,
    step_name=f"trebl_experiment_{name_only}", 
    seq_file=file_path,
    bc_objects=bc_objects,
    umi_object=UMI,
    reverse_complement=True,
    design_file_path=None
)
umi_mapper.create_map()

# Refinement

figures_path = os.path.join(yeast_pool_C_umi_output_path, "figures", name_only)
os.makedirs(figures_path, exist_ok=True)

refiner = map_refiner.MapRefiner(
    db_path=db_path,
    bc_objects=bc_objects,
    column_pairs=[],
    map_order=['quality'],
    step_name=f"trebl_experiment_{name_only}", 
    descriptor="",
    output_figures_path=figures_path,
    reads_threshold=0,
    umi_object=UMI
)
refiner.refine_map_from_db(should_check_exists=True)
refiner.plot_loss()

# Deduplication
deduplicator = umi_deduplicate.UMIDeduplicator(
    db_path=db_path,
    bc_objects=bc_objects,
    step_name=f"trebl_experiment_{name_only}", 
    descriptor="",
    step1_map_name=None,
    fastq_path=file_path,
    output_path=umi_path, 
    refined_map_suffix='quality'
)
deduplicator.run_simple_deduplication()
deduplicator.save_simple_deduplication()
