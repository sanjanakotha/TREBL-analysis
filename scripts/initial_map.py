import pandas as pd
from preprocess import *
import importlib
import error_correct
import duckdb
from pathlib import Path


class InitialMapper:
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

    def __init__(self, db_path, seq_file, bc_objects, step_name, reverse_complement, design_file_path=None, umi_length=0, test_n_reads = None):
        self.con = duckdb.connect(db_path)
        # ensure seq_file is always a list
        self.seq_files = seq_file if isinstance(seq_file, list) else [seq_file]
        self.design_file_path = str(Path(design_file_path).resolve()) if design_file_path else None
        self.bc_objects = bc_objects
        self.reverse_complement = reverse_complement
        self.seq_df = None
        self.mapped_df = None
        self.umi_length = umi_length
        self.step_name = step_name

        self.cols = [bc_object.name for bc_object in bc_objects]
        cols_str = "_".join(self.cols)
        self.table_prefix = f"{step_name}_{cols_str}_"
        self.test_n_reads = test_n_reads

    @time_it
    def read_fastq(self):
        """Read sequences from FASTQ/TXT files and create initial seq table"""
        con = self.con
        if isinstance(self.seq_files, str):
            self.seq_files = [self.seq_files]

        print(f"Reading {len(self.seq_files)} FASTQ/TXT file(s)...")
        seq_queries = []
        for f in self.seq_files:
            f_path = str(Path(f).resolve())
            seq_queries.append(f"""
                SELECT column0 AS sequence
                FROM (
                    SELECT *, ROW_NUMBER() OVER () AS rn
                    FROM read_csv('{f_path}', header=False)
                ) t
                WHERE rn % 4 = 2
            """)

        combined_query = " UNION ALL ".join(seq_queries)

        if self.test_n_reads:
            # Limit to test_n_reads
            combined_query = f"SELECT * FROM ({combined_query}) LIMIT {self.test_n_reads}"

        con.execute(f"CREATE OR REPLACE TABLE seq AS {combined_query};")
        
    @time_it
    def reverse_complement_fastq(self):
        """Replace sequence column with its reverse complement"""
        if not self.reverse_complement:
            return
        print("Reverse complement of sequences...")
        con = self.con
    
        # Create a new table with only the reverse complemented sequence
        con.execute("""
            CREATE TABLE seq_rc AS
            SELECT reverse(translate(sequence, 'ACGTacgt', 'TGCAtgca')) AS sequence
            FROM seq;
        """)
    
        # Drop the old table and rename
        con.execute("DROP TABLE seq;")
        con.execute("ALTER TABLE seq_rc RENAME TO seq;")

    @time_it
    def extract_barcodes(self):
        """Extract barcodes and compute barcode quality"""
        con = self.con
        print(f"Extracting {len(self.bc_objects)} barcodes...")
    
        for bc in self.bc_objects:
            # Build regex for barcode extraction
            if bc.post == "" or bc.preceder == "":
                regex = f"{bc.preceder}(.{{1,{bc.length}}}){bc.post}"
            else:
                regex = f"{bc.preceder}(.*){bc.post}"

            print(regex)
    
            # Create barcode column and fill with regex match (or empty string if none)
            con.execute(f"""
                ALTER TABLE seq ADD COLUMN {bc.name} TEXT;
                UPDATE seq
                SET {bc.name} = coalesce(regexp_extract(sequence::VARCHAR, '{regex}', 1), '');
            """)
    
            # Add quality flag: True if extracted length == expected length
            qual_col = f"{bc.name}_qual"
            con.execute(f"ALTER TABLE seq ADD COLUMN {qual_col} BOOLEAN;")
            con.execute(f"""
                UPDATE seq
                SET {qual_col} = LENGTH({bc.name}) = {bc.length};
            """)
    @time_it
    def extract_umi(self):
        """Extract UMI sequences from the read depending on reverse_complement"""
        print(f"Extracting UMIs ({self.umi_length} bases)...")
        if self.umi_length <= 0:
            return
        con = self.con
    
        con.execute("ALTER TABLE seq ADD COLUMN UMI TEXT;")
    
        if self.reverse_complement:
            # UMI is last self.umi_length bases
            con.execute(f"""
                UPDATE seq SET UMI = right(sequence, {self.umi_length});
            """)
        else:
            # UMI is first self.umi_length bases
            con.execute(f"""
                UPDATE seq SET UMI = left(sequence, {self.umi_length});
            """)
            
    @time_it
    def merge_design(self):
        """Merge with design file if provided, or create default Designed column"""
        con = self.con
        if self.design_file_path:
            print("Merging with design file...")
            con.execute(f"""
                CREATE OR REPLACE TABLE design AS
                SELECT CAST(column0 AS VARCHAR) AS AD
                FROM read_csv_auto('{self.design_file_path}', header=False)
            """)
            con.execute(f"""
                CREATE OR REPLACE TABLE {self.table_prefix}initial AS
                SELECT s.*, CASE WHEN d.AD IS NOT NULL THEN 1 ELSE 0 END AS Designed
                FROM seq s
                LEFT JOIN design d USING(AD);
            """)
        else:
            con.execute(f"""
                CREATE OR REPLACE TABLE {self.table_prefix}initial AS
                SELECT *, 1 AS Designed
                FROM seq;
            """)
        con.execute("DROP TABLE IF EXISTS seq;")

        # Drop the sequence column from the new table
        con.execute(f"ALTER TABLE {self.table_prefix}initial DROP COLUMN sequence;")

    def create_map(self):
        """Master function to run all steps save initial map in duckdb database"""
        self.read_fastq()
        self.reverse_complement_fastq()
        self.extract_barcodes()
        self.extract_umi()
        self.merge_design()
        print("Mapping complete.")

    def preview_map(self):
        print(f"{self.table_prefix}initial")

        # Count total rows
        total_rows = self.con.execute(f"SELECT COUNT(*) FROM {self.table_prefix}initial").fetchone()[0]
        print(f"Total rows: {total_rows}")
        
        return self.con.execute(f"SELECT * FROM {self.table_prefix}initial LIMIT 5;").df()

    
    # def apply_whitelist(self, output_dir):
    #     output_dir = Path(output_dir)
    #     output_dir.mkdir(exist_ok=True)
        
    #     # Run umi_tools whitelist in parallel for all barcodes
    #     error_correct.run_whitelist_parallel(
    #         barcodes=self.bc_objects,
    #         input_fastq=self.seq_files,
    #         output_dir=output_dir,
    #         reverse_complement=self.reverse_complement
    #     )
    
    #     error_correct.plot_all_whitelists(
    #         barcodes=self.bc_objects,
    #         output_dir=output_dir,
    #         reverse_complement=self.reverse_complement
    #     ) 
        
    #     # Process each barcode individually
    #     for bc in self.bc_objects:
    #         # Load whitelist mapping
    #         mapping_df = error_correct.convert_txt_to_whitelist_mapping_df(
    #             bc, output_dir, reverse_complement=self.reverse_complement
    #         )
    
    #         # Push mapping to DuckDB
    #         self.con.execute("CREATE OR REPLACE TABLE barcode_map AS SELECT * FROM mapping_df;")
    
    #         # First, filter initial table to rows where barcode exists in whitelist
    #         self.con.execute(f"""
    #             CREATE OR REPLACE TABLE {self.table_prefix}filtered AS
    #             SELECT m.*
    #             FROM {self.table_prefix}initial AS m
    #             JOIN barcode_map AS b
    #             ON m.{bc.name} = b.original;
    #         """)
    
    #         # Then, update barcode column to canonical values
    #         self.con.execute(f"""
    #             UPDATE {self.table_prefix}filtered AS m
    #             SET {bc.name} = b.canonical
    #             FROM barcode_map AS b
    #             WHERE m.{bc.name} = b.original;
    #         """)
    
    #     # Replace initial table with filtered one
    #     self.con.execute(f"DROP TABLE {self.table_prefix}initial;")
    #     self.con.execute(f"ALTER TABLE {self.table_prefix}filtered RENAME TO {self.table_prefix}initial;")
    
    