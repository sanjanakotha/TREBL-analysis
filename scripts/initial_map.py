import dask.dataframe as dd
import pandas as pd
import subprocess
from Bio.Seq import Seq
from dask.diagnostics import ProgressBar
import os
from finder import add_barcode, add_multiple_barcodes
from preprocess import *

class BarcodeMapper:
    """
    A class to extract and map barcodes from DNA sequences in FASTQ or TXT files.

    Args:
        seq_file (str or list[str]): Path(s) to the input FASTQ or sequence file(s).
        design_file_path (str): Path to the CSV file containing designed sequences.
        bc_names (list of str): Names of barcodes to extract. 
        preceders (list of str): Preceding sequences used for barcode extraction.
        posts (list of str): Following sequences used for barcode extraction.
        lengths (list of int): Lengths of each barcode to extract.
        reverse_complement (bool): Whether to reverse complement sequences before mapping.
        txt_path (str): Path for temporary shortened FASTQ sequences as text.

    Example:
        >>> mapper = BarcodeMapper(
        ...     seq_file="reads.fastq",
        ...     design_file_path="design.csv",
        ...     bc_names=["AD", "AD_BC", "RPTR_BC"],
        ...     preceders=["GGCTAGC", "CGCGCC", "CTCGAG"],
        ...     posts=["", "", ""],
        ...     lengths=[120, 11, 14],
        ...     reverse_complement=True,
        ...     txt_path = "reads_shortened.txt"
        ... )
        >>> mapped_df = mapper.create_map()
    """

    def __init__(self, seq_file, bc_names, preceders, posts, lengths, reverse_complement, txt_path="reads_shortened.txt", design_file_path=None):
        self.seq_file = seq_file
        self.design_file_path = design_file_path  # optional
        self.bc_names = bc_names
        self.preceders = preceders
        self.posts = posts
        self.lengths = lengths
        self.reverse_complement = reverse_complement
        self.seq_df = None
        self.mapped_df = None
        self.txt_path = txt_path

        n = len(bc_names)
        if not (len(preceders) == len(posts) == len(lengths) == n):
            raise ValueError(
                f"All barcode-related lists (bc_names, preceders, posts, and lengths) must have the same length. "
                f"\nGot lengths: bc_names={len(bc_names)}, preceders={len(preceders)}, "
                f"posts={len(posts)}, lengths={len(lengths)}"
            )



    def check_designed(self, map_df):
        """
        Marks sequences as 'Designed' based on a provided design file (optional).
        If no design file is provided, returns the input DataFrame unchanged.

        Args:
            map_df (dask.DataFrame): DataFrame with an 'AD' column of sequences to check.

        Returns:
            dask.DataFrame: DataFrame with an additional 'Designed' column (1 if in design file, else 0).

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> df_checked = mapper.check_designed(mapper.seq_df)
        """
        if self.design_file_path is None:
            map_df["Designed"] = 0
            return map_df
            
        design_file = dd.read_csv(self.design_file_path, header=None, names=["AD"], dtype=str)
        design_file["Designed"] = 1

        map_df["AD"] = map_df["AD"].str.strip()
        design_file["AD"] = design_file["AD"].str.strip()

        merged = dd.merge(map_df, design_file, on="AD", how="left")
        merged["Designed"] = merged["Designed"].fillna(0)
        
        return merged

    def create_map(self):
        """
        Processes the sequence file, extracts barcodes, and maps them to designed sequences.

        Returns:
            dask.DataFrame: DataFrame with barcodes mapped to design file.

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> mapped_df = mapper.create_map()
        """
        dfs = []

        # 1. Load seq files as dask dataframe, translate =and use txt intermediate if needed
        all_seq_df = load_and_shorten_files(self.seq_file, self.reverse_complement, self.txt_path)
    
        # 2. Extract all barcodes
        all_seq_df = add_multiple_barcodes(self.bc_names, self.preceders, self.posts, self.lengths, all_seq_df)
    
        # 3. Apply design check
        self.mapped_df = self.check_designed(all_seq_df).drop(columns="sequence")
    
        return self.mapped_df

    

