import dask.dataframe as dd
import pandas as pd
import subprocess
from Bio.Seq import Seq
from dask.diagnostics import ProgressBar
import os
from dask import delayed, compute
import pathlib
from itertools import islice
from tqdm import tqdm
import time
from functools import wraps

def time_it(func):
    """Decorator to print how long a method takes to run."""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        start = time.time()
        result = func(self, *args, **kwargs)
        elapsed = time.time() - start
        if elapsed < 60:
            print(f"Done in {elapsed:.2f} seconds.\n")
        else:
            print(f"Done in {elapsed/60:.2f} minutes.\n")
        return result
    return wrapper


def load_and_shorten_files(seq_files, reverse_complement):
    """
    Load sequence files (FASTQ, gzipped FASTQ, or TXT) as Dask DataFrame.
    For FASTQ files, keeps only every 4th line starting from line 1 (the sequence line).
    
    Args:
        seq_files (list of str or str): List of file paths.
        reverse_complement (bool): Whether to apply reverse complement.
    
    Returns:
        dask.DataFrame: Loaded sequences as a Dask DataFrame with column 'sequence'.
    """
    # Normalize input
    if isinstance(seq_files, str):
        seq_files = [seq_files]

    dfs = []
    for f in seq_files:
        ext = pathlib.Path(f).suffix.lower()
        
        if ext in [".fastq", ".fq", ".gz"]:
            # Read all lines
            df = dd.read_csv(f, header=None, names=["raw"], dtype=str)

            # Keep only rows 1, 5, 9, ... → sequence lines
            df = df.map_partitions(lambda d: d.iloc[1::4])
            df = df.rename(columns={"raw": "sequence"})
        else:
            # Assume plain TXT file with one sequence per line
            df = dd.read_csv(f, header=None, names=["sequence"], dtype=str)

        dfs.append(df)

    # Concatenate all inputs
    all_seq_df = dd.concat(dfs)

    # Optionally apply reverse complement
    if reverse_complement:
        all_seq_df["sequence"] = all_seq_df["sequence"].map_partitions(
            reverse_complement_series,
            meta=("sequence", str)
        )

    return all_seq_df

def reverse_complement_series(s: pd.Series) -> pd.Series:
    """
    Returns the reverse complement of a pandas Series of DNA sequences.

    Args:
        s (pd.Series): Series of DNA sequences.

    Returns:
        pd.Series: Series of reverse complemented sequences.

    Example:
        >>> s = pd.Series(["ATGC", "GATTACA"])
        >>> BarcodeMapper.reverse_complement_series(s)
        0       GCAT
        1    TGTAATC
        dtype: object
    """
    def safe_rc(seq):
        try:
            return str(Seq(str(seq)).reverse_complement())
        except Exception:
            print(seq)
            return ""
    return s.map(safe_rc)

def reverse_complement(seq):
    return str(Seq(str(seq)).reverse_complement())

    
from tqdm import tqdm
import gzip
from itertools import islice

def shorten_seq_file(infile, outfile, chunk_size=4*100000):
    """
    Fastest extraction of sequences (2nd line of each FASTQ record) to a text file.
    Works with plain FASTQ or gzipped FASTQ.
    
    Args:
        infile (str): Input FASTQ or FASTQ.GZ file.
        outfile (str): Output text file (1 sequence per line).
        chunk_size (int): Number of lines to read at once (multiple of 4).
    """
    print(f"Shortening {infile} to {outfile}...\n")
    opener = gzip.open if infile.endswith(".gz") else open

    with opener(infile, "rt") as fin, open(outfile, "w") as fout:
        pbar = tqdm(desc="Processing chunks", unit="chunk")
        while True:
            lines = list(islice(fin, chunk_size))
            if not lines:
                break
            fout.writelines(lines[1::4])  # Write every 2nd line of 4-line record
            pbar.update(1)
        pbar.close()

                
def save_parquet(df, path):
    """
    Saves the mapped barcode DataFrame to a Parquet file with a progress bar.

    Args:
        path (str): Path to save the Parquet file.

    Returns:
        None

    Example:
        >>> mapper = BarcodeMapper(...initialize...)
        >>> mapper.create_map()
        >>> mapper.save_parquet("mapped_barcodes.parquet")
    """
    if os.path.exists(path):
        print(f"Warning: The path '{path}' already exists and will be added to.\n")

    with ProgressBar():
        df.to_parquet(
            path,
            engine='pyarrow',
            write_index=False
        )