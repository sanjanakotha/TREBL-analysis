print("PYTHON STARTED", flush=True)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob
from tqdm import tqdm
import os

sys.path.insert(0, "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/scripts")
import initial_map
import map_refiner
import complexity
import finder
import preprocess
import plotting
import umi_deduplicate

print("Done with imports")

# Input and output paths
yeast_pool_C_umi_RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/*")
yeast_pool_C_umi_output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_umi" 
print(f"Output path: {yeast_pool_C_umi_output_path}")
print(f"Found {len(yeast_pool_C_umi_RT_seq_files)} RT files")
db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_final.db"

# Initialize results containers
complex_RT_results = []
simple_RT_results = []

# Define barcode objects
EC_RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14)
RT_UMI = finder.Barcode(name = "UMI",
                        preceder = "TGTCAC",
                       post = "",
                       length = 12)

for file_path in yeast_pool_C_umi_RT_seq_files:
    
    # Get the file naeme to use for database
    base_name = os.path.basename(file_path)
    name_only = base_name.split('.')[0]
    print(name_only) 

    # Get the file naeme to use for output
    umi_path = os.path.join(yeast_pool_C_umi_output_path, f"trebl_experiment_yeast_pool_C_umi_{name_only}")
    print(umi_path)

    # Extract UMIs and barcodes from reRTs
    umi_mapper = initial_map.InitialMapper(db_path = db_path,
                                       step_name = f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
                                       seq_file = file_path,
                                       bc_objects = [EC_RPTR_BC],
                                        umi_object=RT_UMI,
                                       reverse_complement = True)
    umi_mapper.create_map()

    # Only keep barcodes of correct length
    refiner = map_refiner.MapRefiner(db_path = db_path,
                                        bc_objects=[EC_RPTR_BC],
                                        column_pairs = [],
                                        map_order = ['quality'],
                                        step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
                                        descriptor = "",
                                        should_check_exists = False,
                                        design_check = False,
                                        output_figures_path='/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_umi/figures',
                                        plot_histograms=True,
                                        reads_threshold = 0,
                                        umi_object = RT_UMI)
    refiner.refine_map_from_db()
    refiner.plot_loss()

    # Run both deduplications
    deduplicator = umi_deduplicate.UMIDeduplicator(db_path = db_path,
                                                        bc_objects = [EC_RPTR_BC],
                                                        step_name = f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
                                                        descriptor = "",
                                                        step1_map_name = None,
                                                        fastq_path = file_path,
                                                        output_path = umi_path, 
                                                       refined_map_suffix = 'quality')

    deduplicator.run_both_deduplications()

   