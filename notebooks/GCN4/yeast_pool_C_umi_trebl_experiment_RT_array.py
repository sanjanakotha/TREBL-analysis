# %%
print("PYTHON STARTED", flush=True)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

from pathlib import Path

# ---------- CLI ARG ----------
# argv[1] = RT fastq file
file_path = sys.argv[1]

# ---------- PATHS ----------
TREBL_ROOT = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL"
sys.path.append(TREBL_ROOT)

from trebl_tools import (
    initial_map,
    map_refiner,
    finder,
    umi_deduplicate,
)

# %%
output_root = f"{TREBL_ROOT}/output/GCN4/yeast_pool_C_umi"
duckdb_root = f"{TREBL_ROOT}/duckdb"

os.makedirs(output_root, exist_ok=True)
os.makedirs(duckdb_root, exist_ok=True)

# ---------- FILE NAMES ----------
base_name = os.path.basename(file_path)
name_only = base_name.split(".")[0]
print(f"Processing: {name_only}", flush=True)

# One DB **per RT file**
db_path = os.path.join(duckdb_root, f"GCN4_{name_only}.db")
print(f"DB path: {db_path}", flush=True)

umi_output_path = os.path.join(
    output_root, f"trebl_experiment_yeast_pool_C_umi_{name_only}"
)
print(f"Output path: {umi_output_path}", flush=True)

# ---------- BARCODE OBJECTS ----------
EC_RPTR_BC = finder.Barcode(
    name="RPTR_BC",
    preceder="CTCGAG",
    post="",
    length=14,
)

RT_UMI = finder.Barcode(
    name="UMI",
    preceder="TGTCAC",
    post="",
    length=12,
)

# ---------- INITIAL MAP ----------
umi_mapper = initial_map.InitialMapper(
    db_path=db_path,
    step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}",
    seq_file=file_path,
    bc_objects=[EC_RPTR_BC],
    umi_object=RT_UMI,
    reverse_complement=True,
    design_file_path=None,
)
umi_mapper.create_map()

# ---------- REFINE ----------
refiner = map_refiner.MapRefiner(
    db_path=db_path,
    bc_objects=[EC_RPTR_BC],
    column_pairs=[],
    map_order=["quality"],
    step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}",
    descriptor="",
    output_figures_path=os.path.join(output_root, "figures"),
    reads_threshold=0,
    umi_object=RT_UMI,
)
refiner.refine_map_from_db(should_check_exists=True)
refiner.plot_loss()

# ---------- DEDUPLICATION ----------
deduplicator = umi_deduplicate.UMIDeduplicator(
    db_path=db_path,
    bc_objects=[EC_RPTR_BC],
    step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}",
    descriptor="",
    step1_map_name=None,
    fastq_path=file_path,
    output_path=umi_output_path,
    refined_map_suffix="quality",
)

deduplicator.run_simple_deduplication()
deduplicator.save_simple_deduplication()

print("DONE", flush=True)
