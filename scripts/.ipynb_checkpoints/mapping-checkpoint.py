import dask.dataframe as dd
import pandas as pd
import subprocess
from Bio.Seq import Seq
from dask.diagnostics import ProgressBar
import os

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

    def __init__(self, seq_file, design_file_path, bc_names, preceders, posts, lengths, reverse_complement, txt_path = "reads_shortened.txt"):
        self.seq_file = seq_file
        
        # Normalize input into list
        if isinstance(self.seq_file, str):
            self.seq_files = [self.seq_file]
        elif isinstance(self.seq_file, (list, tuple)):
            self.seq_files = list(self.seq_file)
        else:
            raise ValueError("seq_file must be a string path or a list/tuple of paths")
    
        self.design_file_path = design_file_path
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

        if "fastq" in seq_file:
            print("fastq file detected. Will shorten to txt file and save at " + txt_path)

    @staticmethod
    def shorten_seq_file(infile, outfile):
        """
        Extracts sequences from a FASTQ or gzipped FASTQ file into a plain text file
        (1 sequence per line). Only the 2nd line of each FASTQ record is extracted.

        Args:
            infile (str): Path to the input FASTQ file (.fastq, .fq, or .gz).
            outfile (str): Path to save the extracted sequences as a text file.

        Returns:
            None

        Example:
            >>> BarcodeMapper.shorten_seq_file(
            ...     infile="reads.fastq.gz",
            ...     outfile="reads_shortened.txt"
            ... )
        """
        
        import gzip
    
        opener = gzip.open if infile.endswith(".gz") else open
    
        with opener(infile, "rt") as fin, open(outfile, "w") as fout:
            for i, line in enumerate(fin):
                if i % 4 == 1:   # 2nd line of each FASTQ record
                    fout.write(line)
                    
    @staticmethod
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

    @staticmethod
    def add_barcode(seq_df, bc_name, preceder, post, bc_length=120):
        """
        Extracts a barcode from sequences between a preceding and following pattern.

        Args:
            seq_df (pd.DataFrame or dask.DataFrame): DataFrame containing a 'sequence' column.
            bc_name (str): Name of the barcode to create.
            preceder (str): Sequence preceding the barcode.
            post (str): Sequence following the barcode.
            bc_length (int, optional): Maximum length of barcode to extract. Defaults to 120.

        Returns:
            pd.DataFrame or dask.DataFrame: Updated DataFrame with barcode and quality columns.

        Example:
            >>> df = pd.DataFrame({"sequence": ["AAAXXXCCC", "AAAYYYCCC"]})
            >>> BarcodeMapper.add_barcode(df, "AD", "AAA", "CCC", 3)
               sequence   AD  AD_qual
            0  AAAXXXCCC  XXX     True
            1  AAAYYYCCC  YYY     True
        """
        regex = f"{preceder}(.*){post}"
        subseq_series = seq_df['sequence'].str.extract(regex)[0].str.slice(0, bc_length)
        seq_df[bc_name] = subseq_series
        seq_df[bc_name + "_qual"] = seq_df[bc_name].str.len() == bc_length
        seq_df[bc_name + "_qual"] = seq_df[bc_name + "_qual"].fillna(False)
        return seq_df

    def add_multiple_barcodes(self, seq_df):
        """
        Extracts all defined barcodes from sequences in a DataFrame.

        Args:
            seq_df (pd.DataFrame or dask.DataFrame): DataFrame containing sequences.

        Returns:
            pd.DataFrame or dask.DataFrame: Updated DataFrame with all barcodes added.

        Example:
            >>> df = pd.DataFrame({"sequence": ["AAAXXXCCC", "AAAYYYCCC"]})
            >>> mapper = BarcodeMapper(...initialize...)
            >>> mapper.add_multiple_barcodes(df)
        """
        df = seq_df.copy()
        for bc_name, preceder, post, length in zip(self.bc_names, self.preceders, self.posts, self.lengths):
            df = self.add_barcode(df, bc_name, preceder, post, length)
        return df

    def check_designed(self, map_df):
        """
        Marks sequences as 'Designed' based on a provided design file.

        Args:
            map_df (dask.DataFrame): DataFrame with an 'AD' column of sequences to check.

        Returns:
            dask.DataFrame: DataFrame with an additional 'Designed' column (1 if in design file, else 0).

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> df_checked = mapper.check_designed(mapper.seq_df)
        """
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
        processed_dfs = []
    
        for i, f in enumerate(self.seq_files):
            ext = os.path.splitext(f)[1].lower()
            input_file = f
        
            if ext in [".fastq", ".fq", ".gz"]:   # also handle gzipped
                out_path = self.txt_path if len(self.seq_files) == 1 else f"{self.txt_path}_{i}.txt"
                self.shorten_seq_file(f, out_path)   # <-- pass just this file
                input_file = out_path
    
            # Load sequences into DataFrame
            seq_df = dd.read_csv(input_file, header=None, names=["sequence"], dtype=str)
    
            # Reverse complement if requested
            if self.reverse_complement:
                seq_df['sequence'] = seq_df['sequence'].map_partitions(
                    self.reverse_complement_series,
                    meta=('sequence', str)
                )
    
            # Extract barcodes + check against design
            seq_df = self.add_multiple_barcodes(seq_df)
            mapped_df = self.check_designed(seq_df)
            mapped_df = mapped_df.drop(columns="sequence")
    
            processed_dfs.append(mapped_df)
    
        # Concatenate all results
        self.mapped_df = dd.concat(processed_dfs)
    
        return self.mapped_df

    def save_parquet(self, path):
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
            self.mapped_df.to_parquet(
                path,
                engine='pyarrow',
                write_index=False
            )

