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

EC_AD = finder.Barcode(name = "AD",
                       preceder = "GGCTAGC",
                       post = "TGACTAG",
                       length = 120)

EC_AD_BC = finder.Barcode(name = "AD_BC",
                       preceder = "CGCGCC",
                       post = "GGGCCC",
                       length = 11)

EC_RPTR_BC = finder.Barcode(name = "RPTR_BC",
                       preceder = "CTCGAG",
                       post = "GGCCGC",
                       length = 14)


step1_mapper = initial_map.InitialMapper(db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/GCN4_with_post.db",
                                       step_name = "step1", #Should include time point, replicate (don't need to specify BC type(s))
                                       seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A10_S1.fastq.gz.assembled.fastq",
                                       design_file_path =  "/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/A10_sequencing/v2/current/a10_designfile.csv",
                                       bc_objects = [EC_AD, EC_AD_BC, EC_RPTR_BC],
                                       reverse_complement = True)
step1_mapper.create_map()

