import dask.dataframe as dd
import pandas as pd
import subprocess
from Bio.Seq import Seq
from dask.diagnostics import ProgressBar


class BarcodeMapper:
    """
    A class to extract and map barcodes from DNA sequences in FASTQ or CSV files.

    Args:
        seq_file (str): Path to the input FASTQ or sequence file.
        design_file_path (str): Path to the CSV file containing designed sequences.
        bc_names (list of str): Names of barcodes to extract. 
        preceders (list of str): Preceding sequences used for barcode extraction.
        posts (list of str): Following sequences used for barcode extraction.
        lengths (list of int): Lengths of each barcode to extract.
        reverse_complement (bool): Whether to reverse complement sequences before mapping.
    Example:
        >>> mapper = BarcodeMapper(
        ...     seq_file="reads.fastq",
        ...     design_file_path="design.csv",
        ...     bc_names=["AD", "AD_BC", "RPTR_BC"],
        ...     preceders=["GGCTAGC", "CGCGCC", "CTCGAG"],
        ...     posts=["", "", ""],
        ...     lengths=[120, 11, 14]
        ...     reverse_complement=True
        ... )
        >>> mapped_df = mapper.create_map()
    """

    def __init__(self, seq_file, design_file_path, bc_names, preceders, posts, lengths, reverse_complement):
        self.seq_file = seq_file
        self.design_file_path = design_file_path
        self.bc_names = bc_names
        self.preceders = preceders
        self.posts = posts
        self.lengths = lengths
        self.reverse_complement = reverse_complement
        self.seq_df = None
        self.mapped_df = None

        n = len(bc_names)
        if not (len(preceders) == len(posts) == len(lengths) == n):
            raise ValueError(
                f"All barcode-related lists (bc_names, preceders, posts, and lengths) must have the same length. "
                f"\nGot lengths: bc_names={len(bc_names)}, preceders={len(preceders)}, "
                f"posts={len(posts)}, lengths={len(lengths)}"
            )

    def shorten_seq_file(self, output_file):
        """
        Extracts the DNA sequences (2nd line of every FASTQ record) using awk.

        Args:
            output_file (str): Path to save the extracted sequences.

        Returns:
            None

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> mapper.shorten_seq_file("short_reads.txt")
        """
        cmd = f"awk 'NR%4 == 2' {self.seq_file} > {output_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Extraction complete. DNA sequences written to {output_file}")
        self.seq_file = output_file

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

    def create_map(self, should_shorten_seq_file=False):
        """
        Processes the sequence file, extracts barcodes, and maps them to designed sequences.

        Args:
            should_shorten_seq_file (str or bool, optional): If provided as a path, the
                sequence file is first shortened and saved to this path. Defaults to False.

        Returns:
            dask.DataFrame: DataFrame with barcodes mapped to design file.

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> mapped_df = mapper.create_map()
        """
        if should_shorten_seq_file:
            self.shorten_seq_file(should_shorten_seq_file)
            self.seq_file = should_shorten_seq_file
        
        self.seq_df = dd.read_csv(self.seq_file, header=None, names=["sequence"], dtype=str)

        if self.reverse_complement:
            self.seq_df['sequence'] = self.seq_df['sequence'].map_partitions(
                self.reverse_complement_series,
                meta=('sequence', str)
            )

        self.seq_df = self.add_multiple_barcodes(self.seq_df)        
        self.mapped_df = self.check_designed(self.seq_df)
        self.mapped_df = self.mapped_df.drop(columns="sequence")
        
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
        with ProgressBar():
            self.mapped_df.to_parquet(
                path,
                engine='pyarrow',
                write_index=False
            )

