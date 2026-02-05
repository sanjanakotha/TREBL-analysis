print("PYTHON STARTED", flush=True)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob
from tqdm import tqdm
import os

sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/")
from scripts import initial_map, map_refiner, complexity, finder, preprocess, error_correct, plotting, umi_deduplicate

print("Done with imports", flush=True)

# Input and output paths
yeast_pool_C_umi_RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/*fastq*")
#yeast_pool_C_umi_RT_seq_files = ["global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/RPTR_3_30_S19_R1_001_fastp.fastq.gz"]

yeast_pool_C_umi_output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_umi" 
print(f"Output path: {yeast_pool_C_umi_output_path}", flush=True)
print(f"Found {len(yeast_pool_C_umi_RT_seq_files)} RT files", flush=True)
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
    print(name_only, flush=True)

    # Get the file naeme to use for output
    umi_path = os.path.join(yeast_pool_C_umi_output_path, f"trebl_experiment_yeast_pool_C_umi_{name_only}")
    print(umi_path, flush=True)

    # Extract UMIs and barcodes from reRTs
    umi_mapper = initial_map.InitialMapper(db_path = db_path,
                                       step_name = f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
                                       seq_file = file_path,
                                       bc_objects = [EC_RPTR_BC],
                                        umi_object=RT_UMI,
                                       reverse_complement = True,
                                          design_file_path = None)
    umi_mapper.create_map()

    # Only keep barcodes of correct length
    refiner = map_refiner.MapRefiner(db_path = db_path,
                                        bc_objects=[EC_RPTR_BC],
                                        column_pairs = [],
                                        map_order = ['quality'],
                                         # map_order = ['grouped'],
                                        step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}", 
                                        descriptor = "",
                                        output_figures_path='/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_umi/figures',
                                        reads_threshold = 0,
                                        umi_object = RT_UMI)
    refiner.refine_map_from_db(should_check_exists=True)
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

    #deduplicator.run_both_deduplications()
    deduplicator.run_simple_deduplication()
    deduplicator.save_simple_deduplication()


