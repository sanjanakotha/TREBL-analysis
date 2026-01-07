import sys, os
sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/scripts")

import finder
import umi_deduplicate
from tqdm import tqdm  # progress bar

EC_AD_BC = finder.Barcode(name = "AD_BC",
                       preceder = "CGCGCC",
                       post = "GGGCCC",
                       length =11)

full_deduplicator = umi_deduplicate.UMIToolsDeduplicator(db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4.db",
                                                    bc_objects = [EC_AD_BC],
                                                    step_name = "umi_1_15",
                                                    descriptor = "",
                                                    step1_map_name = "",
                                                    fastq_path = "/global/scratch/projects/fc_mvslab/data/sequencing/20250218_MZCCSCU_MedGenome/MZ/results/assembled/AD_1_15.fastq.gz.assembled.fastq",
                                                    umi_length = 12,
                                                    reverse_complement = True,
                                                    output_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/AD_1_15")
full_deduplicator.run_full_umi_pipeline()