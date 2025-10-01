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

def load_and_shorten_files(seq_files, reverse_complement, txt_path="reads_shortened.txt"):
    """
    Load sequence files (FASTQ, gzipped FASTQ, or TXT) as dask dataframe and shorten FASTQ files if needed.
    
    Args:
        seq_files (list of str): List of file paths.
        reverse_complement (bool): Whether to apply reverse complement.
        txt_path (str): Base path for shortened text files.

    Returns:
        dask.DataFrame: Loaded sequences as a Dask DataFrame.
    """
    # Normalize input
    if isinstance(seq_files, str):
        seq_files = [seq_files]

    @delayed
    def shorten_if_needed(f, out_path):
        ext = pathlib.Path(f).suffix.lower()
        
        try:
            if ext in [".fastq", ".fq", ".gz"]:
                print(f"Shortening {f} -> {out_path}")
                shorten_seq_file(f, out_path)
            else:
                out_path = f  # already a TXT file
            return out_path
        except Exception as e:
            print(f"Error processing file {f}: {e}")
            raise

    # Launch parallel tasks
    delayed_files = [
        shorten_if_needed(f, txt_path if len(seq_files) == 1 else f"{txt_path.rsplit('.', 1)[0]}_{i}.txt"
)
        for i, f in enumerate(seq_files)
    ]

    # Compute shortened files in parallel
    shortened_files = compute(*delayed_files)

    print(shortened_files)

    for path in shortened_files:
        if not os.path.exists(path):
            print(f"Missing: {path}")
        elif os.path.getsize(path) == 0:
            print(f"Empty: {path}")

    # Load all TXT files into a single Dask DataFrame
    all_seq_df = dd.read_csv(shortened_files, header=None, names=["sequence"], dtype=str)

    if reverse_complement:
        all_seq_df['sequence'] = all_seq_df['sequence'].map_partitions(
            reverse_complement_series,
            meta=('sequence', str)
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
    

def shorten_seq_file(infile, outfile):
    """
    Extracts sequences from a FASTQ/FASTQ.GZ file into a plain text file
    (1 sequence per line). Only the 2nd line of each FASTQ record is extracted.
    """
    opener = gzip.open if infile.endswith(".gz") else open
    
    with opener(infile, "rt") as fin, open(outfile, "w") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 1:  # 2nd line of each FASTQ record
                fout.write(line)

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
        print(f"Warning: The path '{path}' already exists and will be added to.")

    with ProgressBar():
        df.to_parquet(
            path,
            engine='pyarrow',
            write_index=False
        )