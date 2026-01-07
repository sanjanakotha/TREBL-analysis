import sys
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
import dask.dataframe as dd
from tqdm import tqdm  # progress bar

sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/scripts")
import initial_map
import map_refiner
import complexity
import finder
import preprocess
import plotting
import umi_deduplicate

# With BC Post
print("First, with post")

AD = finder.Barcode(name="AD",
                             preceder="AGGAGCA",
                             post="TGATAAG",
                             length=186)

AD_BC = finder.Barcode(name="AD_BC",
                             preceder="GGCCTC",
                             post="GGGCCC",
                             length=20)
RPTR_BC = finder.Barcode(name="RPTR_BC",
                             preceder="CTCGAG",
                             post="GGCCGC",
                             length=20)

step1_mapper = initial_map.InitialMapper(db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/Plasm.db",
                                       step_name = "step1", 
                                       seq_file =  "/global/scratch/projects/fc_mvslab/data/sequencing/czb_new_sept2025/Marissa/results/Plasm_step1_S1.fastq.gz.assembled.fastq",
                                       design_file_path =  None,
                                       bc_objects = [AD, AD_BC, RPTR_BC],
                                       reverse_complement = False)

step1_mapper.create_map()


# Without BC Post
print("Next, without post")

AD = finder.Barcode(name="AD",
                             preceder="AGGAGCA",
                             post="",
                             length=186)

AD_BC = finder.Barcode(name="AD_BC",
                             preceder="GGCCTC",
                             post="",
                             length=20)

RPTR_BC = finder.Barcode(name="RPTR_BC",
                             preceder="CTCGAG",
                             post="",
                             length=20)

step1_mapper = initial_map.InitialMapper(db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/Plasm_no_post.db",
                                       step_name = "step1", 
                                       seq_file =  "/global/scratch/projects/fc_mvslab/data/sequencing/czb_new_sept2025/Marissa/results/Plasm_step1_S1.fastq.gz.assembled.fastq",
                                       design_file_path =  None,
                                       bc_objects = [AD, AD_BC, RPTR_BC],
                                       reverse_complement = False)

step1_mapper.create_map()
