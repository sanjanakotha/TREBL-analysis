import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob
sys.path.append("../../scripts")  # allow imports from local scripts directory
import initial_map
import map_refiner
import complexity
import finder
import preprocess
import plotting
import os
from tqdm import tqdm  # progress bar

print("Done with imports")

# Define barcode objects for AD and reporter reads
EC_AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120)
EC_AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11)
EC_RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14)

# Path to DuckDB database
db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_final.db"
print(f"Using database: {db_path}")

# Input FASTQ files for AD and RPTR reads
yeast_pool_C_PolyT_AD_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Sep2024/ciber2_iii_MZ001/EC_Ciber2_iii_Gcn4/results/assembled/AD/*")
yeast_pool_C_PolyT_RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Sep2024/ciber2_iii_MZ001/EC_Ciber2_iii_Gcn4/results/assembled/RPTR/*")
print(f"Found {len(yeast_pool_C_PolyT_AD_seq_files)} AD files")
print(f"Found {len(yeast_pool_C_PolyT_RT_seq_files)} RPTR files")

# Output directory for processed results
yeast_pool_C_PolyT_output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_PolyT" 
os.makedirs(yeast_pool_C_PolyT_output_path, exist_ok=True)
print(f"Output path: {yeast_pool_C_PolyT_output_path}")

# Group AD barcodes for mapping
AD_objects = [EC_AD, EC_AD_BC]

# Process AD sequencing files
for file_path in yeast_pool_C_PolyT_AD_seq_files:
    base_name = os.path.basename(file_path)  # extract filename
    name_only = base_name.split('.')[0].replace("-", "_")  # sanitize step name
    print(f"Processing AD file: {name_only}")

    file_output_path = os.path.join(yeast_pool_C_PolyT_output_path, f"trebl_experiment_yeast_pool_C_PolyT_{name_only}")
    print(f"Output prefix: {file_output_path}")

    # Map AD barcodes and UMIs from reads
    bc_mapper = initial_map.InitialMapper(
        db_path=db_path,
        step_name=f"trebl_experiment_yeast_pool_C_PolyT_{name_only}",
        seq_file=file_path,
        design_file_path="/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/A10_sequencing/v2/current/a10_designfile.csv",
        bc_objects=AD_objects,
        reverse_complement=True
    )
    print("Creating AD barcode map")
    bc_mapper.create_map()

    # Refine AD barcode map with design and read-count filtering
    refiner = map_refiner.MapRefiner(
        db_path=db_path,
        bc_objects=AD_objects,
        column_pairs=[],
        map_order=['quality', 'designed', 'grouped', 'thresholded'],
        step_name=f"trebl_experiment_yeast_pool_C_PolyT_{name_only}",
        descriptor="",
        should_check_exists=False,
        design_check=True,
        reads_threshold=20
    )
    print("Refining AD barcode map")
    refiner.refine_map_from_db()

# Group RPTR barcodes for mapping
RPTR_objects = [EC_RPTR_BC]

# Process RPTR sequencing files
for file_path in yeast_pool_C_PolyT_RT_seq_files:
    base_name = os.path.basename(file_path)  # extract filename
    name_only = base_name.split('.')[0].replace("-", "_")  # sanitize step name
    print(f"Processing RPTR file: {name_only}")

    file_output_path = os.path.join(yeast_pool_C_PolyT_output_path, f"trebl_experiment_yeast_pool_C_PolyT_{name_only}")
    print(f"Output prefix: {file_output_path}")

    # Map RPTR barcodes from reads
    bc_mapper = initial_map.InitialMapper(
        db_path=db_path,
        step_name=f"trebl_experiment_yeast_pool_C_PolyT_{name_only}",
        seq_file=file_path,
        design_file_path=None,
        bc_objects=RPTR_objects,
        reverse_complement=True
    )
    print("Creating RPTR barcode map")
    umi_map = bc_mapper.create_map()

    # Refine RPTR barcode map with read-count filtering
    refiner = map_refiner.MapRefiner(
        db_path=db_path,
        bc_objects=RPTR_objects,
        column_pairs=[],
        map_order=['quality', 'grouped', 'thresholded'],
        step_name=f"trebl_experiment_yeast_pool_C_PolyT_{name_only}",
        descriptor="",
        should_check_exists=False,
        design_check=False,
        reads_threshold=20
    )
    print("Refining RPTR barcode map")
    refiner.refine_map_from_db()