import sys
import glob
import os
sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/")
from scripts import finder, pipelines
from tqdm import tqdm
print("imports done")

# Initialize TREBL Pipeline
pipeline_no_err_corr = pipelines.TreblPipeline(
    db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/NKX2-2_TL4B2_no_err_corr_rim.db", # CHANGE TO DUCKDB LOCATION
    design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/DNA_Tiles_nkx2_2.txt", 
    error_correction = False,  
    output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/NKX2-2/TL4B2/rim/no_err_corr")  # CHANGE TO DESIRED OUTPUT, FOLDER MUST EXIST

# Input sequence files
AD_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Caitlin/TL4B2/rim/a*fastq") # CHANGE TO CORRECT PATHS
RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Caitlin/TL4B2/rim/r*fastq") # CHANGE TO CORRECT PATHS

# Define barcodes
trebl_exp_ADBC2 = finder.Barcode(name="ADBC2", preceder="TATGCTAT", post="GGCCGGCCG", length=6)
trebl_exp_HawkBCs = finder.Barcode(name="HawkBCs", preceder="TAGC", post="CTCGAGA", length=9)
trebl_exp_RTBC = finder.Barcode(name="RTBC", preceder="GCCCC", post="GCGG", length=16)

# UMI barcode (we just want last 12 positions so only specifying length here)
UMI = finder.Barcode(name="UMI", preceder="", post="", length=12)

print("initialized pipeline and saved seq files and barcodes")

# Run TREBL experiment analysis
pipeline_no_err_corr.trebl_experiment_analysis(
        AD_seq_files = AD_seq_files,
        AD_bc_objects = [trebl_exp_ADBC2, trebl_exp_HawkBCs],
        RT_seq_files = RT_seq_files,
        RT_bc_objects = [trebl_exp_RTBC],
        reverse_complement = True, # Should you reverse complement sequences?
        AD_umi_object = UMI,       # UMI object for AD reads
        RT_umi_object = UMI,       # UMI object for RT reads
        reads_threshold_AD=0,      # Threshold of 1 because no EC
        reads_threshold_RT=0,       # Threshold of 1 because no EC
    )