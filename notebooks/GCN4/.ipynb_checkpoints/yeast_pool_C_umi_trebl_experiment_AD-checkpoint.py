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
from tqdm import tqdm
import umi_deduplicate

print("Done with imports")

# Input and output paths
yeast_pool_C_umi_AD_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/20250218_MZCCSCU_MedGenome/MZ/results/assembled/*")
yeast_pool_C_umi_output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_umi" 
os.makedirs(yeast_pool_C_umi_output_path, exist_ok=True)
print(f"Output path: {yeast_pool_C_umi_output_path}")
print(f"Found {len(yeast_pool_C_umi_AD_seq_files)} AD files")
db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_final.db"

# Initialize results containers
complex_AD_results = []
simple_AD_results = []

# Define barcode objects
EC_AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120)
EC_AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11)
EC_RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14)
AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)

# AD barcodes for mapping
AD_objects = [EC_AD, EC_AD_BC]

# Process each AD sequencing file
for file_path in yeast_pool_C_umi_AD_seq_files:
    base_name = os.path.basename(file_path)  # filename
    name_only = base_name.split('.')[0]  # step name
    print(f"Processing file: {name_only}")

    umi_path = os.path.join(yeast_pool_C_umi_output_path, f"trebl_experiment_yeast_pool_C_umi_{name_only}")
    print(f"Output prefix: {umi_path}")

    # Map AD barcodes and UMIs
    umi_mapper = initial_map.InitialMapper(
        db_path=db_path,
        step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}",
        seq_file=file_path,
        design_file_path="/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/A10_sequencing/v2/current/a10_designfile.csv",
        bc_objects=AD_objects,
        umi_object=AD_UMI,
        reverse_complement=True
    )
    print("Creating AD barcode map")
    umi_mapper.create_map()

    # Refine AD map
    refiner = map_refiner.MapRefiner(
        db_path=db_path,
        bc_objects=AD_objects,
        column_pairs=[],
        map_order=['quality', 'designed'],
        step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}",
        descriptor="",
        should_check_exists=False,
        design_check=True,
        output_figures_path='/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4/yeast_pool_C_umi/figures',
        plot_histograms=True,
        reads_threshold=0,
        umi_object=AD_UMI
    )
    print("Refining AD barcode map")
    refiner.refine_map_from_db()
    refiner.plot_loss()

    # Run simple and UMI-tools deduplications
    deduplicator = umi_deduplicate.UMIDeduplicator(
        db_path=db_path,
        bc_objects=AD_objects,
        step_name=f"trebl_experiment_yeast_pool_C_umi_{name_only}",
        descriptor="",
        step1_map_name=None,
        fastq_path=file_path,
        output_path=umi_path,
        refined_map_suffix='designed'
    )
    print("Running both deduplications")
    deduplicator.run_both_deduplications()

    # Collect complex results
    one_file_complex_results = pd.read_csv(
        os.path.join(umi_path, f"{name_only}_directional_umi_counts.tsv"), sep="\t"
    )
    one_file_complex_results["name"] = name_only
    complex_AD_results.append(one_file_complex_results)

    # Collect simple results
    one_file_simple_results = pd.read_csv(
        os.path.join(umi_path, f"{name_only}_simple_umi_counts.tsv"), sep="\t"
    )
    one_file_simple_results["name"] = name_only
    simple_AD_results.append(one_file_simple_results)
