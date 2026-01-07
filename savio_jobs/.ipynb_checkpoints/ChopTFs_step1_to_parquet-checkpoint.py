import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
from concurrent.futures import ThreadPoolExecutor, as_completed
import dask.dataframe as dd
import sys

sys.path.append("../scripts")

import initial_map
import map_refiner
import complexity
import finder
import preprocess
import complexity

  
import glob
from tqdm import tqdm  # progress bar


AD = finder.Barcode(name="AD",
                             preceder="GGCTAGC",
                             post="TGACTAG",
                             length=120)

AD_BC = finder.Barcode(name="AD_BC",
                             preceder="CGCGCC",
                             post="GGGCCC",
                             length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC",
                             preceder="CTCGAG",
                             post="GGCCGC",
                             length=14)

# Took about 5 minutes to shorten  "/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A11_S2.fastq.gz.assembled.fastq"
mapper = initial_map.InitialMapper(seq_file =  "../data/A11_S2_reads_shortened.txt",
                        bc_objects=[AD, AD_BC, RPTR_BC],
                      reverse_complement = True,
                                  design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/ChopTFDesign.csv")
mapped_df = mapper.create_map()

# Only had to run once
preprocess.save_parquet(mapped_df, "../data/A11_S2_parquet")