#!/usr/bin/env python3

"""
UMI pipeline converted from Jupyter Notebook for SLURM execution
Author: Caitlin / adapted for batch execution
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
from concurrent.futures import ThreadPoolExecutor, as_completed
import dask.dataframe as dd
import sys
import glob
import os
from Bio.Seq import Seq
from tqdm import tqdm

# ---------------------------------------------------------------------
# Add TREBL codebase to path
# ---------------------------------------------------------------------
sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/")

from scripts import (
    initial_map,
    map_refiner,
    complexity,
    finder,
    preprocess,
    error_correct,
    plotting,
    umi_deduplicate,
)

# ---------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------
def main():

    # -----------------------------------------------------------------
    # Input FASTQ files and output directory
    # -----------------------------------------------------------------
    AD_seq_files = glob.glob(
        "/global/scratch/projects/fc_mvslab/OpenProjects/Caitlin/TL4B2/puro_only/a*fastq"
    )

    RT_seq_files = glob.glob(
        "/global/scratch/projects/fc_mvslab/OpenProjects/Caitlin/TL4B2/puro_only/r*fastq"
    )

    output_path = (
        "/global/scratch/projects/fc_mvslab/OpenProjects/Caitlin/TL4B2/puro_only/results_umi_tools_threshold"
    )

    os.makedirs(output_path, exist_ok=True)

    # -----------------------------------------------------------------
    # Barcode definitions
    # -----------------------------------------------------------------
    step1_ADBC2 = finder.Barcode(
        name="ADBC2",
        preceder="TATGCTAT",
        post="GGCCGGCCG",
        length=6,
    )

    step1_HawkBCs = finder.Barcode(
        name="HawkBCs",
        preceder="TAGC",
        post="CTCGAGA",
        length=9,
    )

    step1_RTBC = finder.Barcode(
        name="RTBC",
        preceder="GCCCC",
        post="GCGG",
        length=16,
    )

    AD_objects = [step1_ADBC2, step1_HawkBCs]
    RTBC_objects = [step1_RTBC]

    # -----------------------------------------------------------------
    # Step 1 map and DuckDB path
    # -----------------------------------------------------------------
    step1_csv = (
        "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/"
        "NKX2-2_whitelist_concat_density_step1_map.csv"
    )

    db_path = os.path.join(
        output_path,
        "NKX2-2_CC_UMI_whitelisted_concat_custom_threshold_redone.db",
    )

    # -----------------------------------------------------------------
    # Containers for results
    # -----------------------------------------------------------------
    complex_AD_results = []
    simple_AD_results = []

    complex_RT_results = []
    simple_RT_results = []

    # =================================================================
    # AD BARCODE FILES
    # =================================================================
    for file_path in AD_seq_files:

        base_name = os.path.basename(file_path)
        name_only = base_name.split(".")[0]
        print(f"\nProcessing AD file: {name_only}")

        umi_path = os.path.join(output_path, f"umi_{name_only}")
        os.makedirs(umi_path, exist_ok=True)

        # --------------------------------------------------------------
        # Initial mapping (UMI + barcodes)
        # --------------------------------------------------------------

        # SK: Commeted out because we don't need to do the initial mapping again
        # umi_mapper = initial_map.InitialMapper(
        #     db_path=db_path,
        #     step_name=f"umi_{name_only}",
        #     seq_file=file_path,
        #     design_file_path=None,
        #     bc_objects=AD_objects,
        #     reverse_complement=True,
        #     umi_object=finder.Barcode(
        #         name="UMI", preceder="", post="", length=12
        #     ),
        # )

        # umi_mapper.create_map()

        # # Print preview instead of display()
        # try:
        #     print(umi_mapper.preview_map().head())
        # except Exception:
        #     print("Preview not available")

        # --------------------------------------------------------------
        # Refinement + error correction
        # --------------------------------------------------------------

         # SK: Changed to use reads threshold of 10 for refining
        refiner = map_refiner.MapRefiner(
            db_path=db_path,
            bc_objects=AD_objects,
            column_pairs=[],
            reads_threshold=10,
            map_order=["quality", "error_corrected"],
            step_name=f"umi_{name_only}",
            descriptor="",
            output_figures_path=umi_path,
            manual_ec_threshold=True,
        )

        refiner.refine_map_from_db()
        refiner.plot_loss()
        refiner.plot_error_correction()

        # --------------------------------------------------------------
        # UMI deduplication
        # --------------------------------------------------------------
        deduplicator = umi_deduplicate.UMIDeduplicator(
            db_path=db_path,
            bc_objects=AD_objects,
            step_name=f"umi_{name_only}",
            descriptor="",
            step1_map_name=None,
            fastq_path=file_path,
            output_path=umi_path,
            refined_map_suffix="error_corrected",
        )

        deduplicator.run_both_deduplications()

        # --------------------------------------------------------------
        # Collect results
        # --------------------------------------------------------------
        complex_df = pd.read_csv(
            os.path.join(
                umi_path, f"{name_only}_directional_umi_counts.tsv"
            ),
            sep="\t",
        )
        complex_df["name"] = name_only
        complex_AD_results.append(complex_df)

        simple_df = pd.read_csv(
            os.path.join(
                umi_path, f"{name_only}_simple_umi_counts.tsv"
            ),
            sep="\t",
        )
        simple_df["name"] = name_only
        simple_AD_results.append(simple_df)

    # =================================================================
    # RT BARCODE FILES
    # =================================================================
    for file_path in RT_seq_files:

        base_name = os.path.basename(file_path)
        name_only = base_name.split(".")[0]
        print(f"\nProcessing RT file: {name_only}")

        umi_path = os.path.join(output_path, f"umi_{name_only}")
        os.makedirs(umi_path, exist_ok=True)

        # --------------------------------------------------------------
        # Initial mapping
        # --------------------------------------------------------------
        umi_mapper = initial_map.InitialMapper(
            db_path=db_path,
            step_name=f"umi_{name_only}",
            seq_file=file_path,
            design_file_path=None,
            bc_objects=RTBC_objects,
            reverse_complement=True,
            umi_object=finder.Barcode(
                name="UMI", preceder="", post="", length=12
            ),
        )

        umi_mapper.create_map()

        # --------------------------------------------------------------
        # Refinement
        # --------------------------------------------------------------
        refiner = map_refiner.MapRefiner(
            db_path=db_path,
            bc_objects=RTBC_objects,
            column_pairs=[],
            reads_threshold=0,
            map_order=["quality", "error_corrected"],
            step_name=f"umi_{name_only}",
            descriptor="",
            output_figures_path=umi_path,
            manual_ec_threshold=False,
        )

        refiner.refine_map_from_db()
        refiner.plot_loss()
        refiner.plot_error_correction()

        # --------------------------------------------------------------
        # UMI deduplication
        # --------------------------------------------------------------
        deduplicator = umi_deduplicate.UMIDeduplicator(
            db_path=db_path,
            bc_objects=RTBC_objects,
            step_name=f"umi_{name_only}",
            descriptor="",
            step1_map_name=None,
            fastq_path=file_path,
            output_path=umi_path,
            refined_map_suffix="error_corrected",
        )

        deduplicator.run_both_deduplications()

        # --------------------------------------------------------------
        # Collect results
        # --------------------------------------------------------------
        complex_df = pd.read_csv(
            os.path.join(
                umi_path, f"{name_only}_directional_umi_counts.tsv"
            ),
            sep="\t",
        )
        complex_df["name"] = name_only
        complex_RT_results.append(complex_df)

        simple_df = pd.read_csv(
            os.path.join(
                umi_path, f"{name_only}_simple_umi_counts.tsv"
            ),
            sep="\t",
        )
        simple_df["name"] = name_only
        simple_RT_results.append(simple_df)

    print("\nPipeline finished successfully.")


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
