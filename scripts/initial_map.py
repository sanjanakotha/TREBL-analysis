from pathlib import Path
import duckdb
from scripts.preprocess import time_it   

class InitialMapper:
    """
    Extract and map barcodes from DNA sequences in FASTQ or TXT files.

    This class reads sequence files, optionally reverse-complements them, 
    extracts barcodes and UMIs, and merges with a design file to produce 
    an initial mapping table in DuckDB.

    Args:
        db_path (str): Path to DuckDB database file.
        seq_file (str or list[str]): Path(s) to input FASTQ or sequence files.
        bc_objects (list): List of barcode objects defining name, preceder, post, and length.
        step_name (str): Prefix for tables in DuckDB.
        reverse_complement (bool): Whether to reverse complement sequences before mapping.
        design_file_path (str, optional): Path to CSV file containing designed sequences. 
        umi_object (object, optional): UMI object defining extraction parameters. Defaults to None.

    Example:
        >>> mapper = InitialMapper(
        ...     db_path="db.duckdb",
        ...     seq_file="reads.fastq",
        ...     bc_objects=[...],
        ...     step_name="step1",
        ...     reverse_complement=True,
        ...     design_file_path="design.csv"
        ... )
        >>> mapper.create_map()
    """

    def __init__(self, 
                 db_path, 
                 seq_file, 
                 bc_objects, 
                 step_name, 
                 reverse_complement, 
                 design_file_path, 
                 umi_object=None):
        """Initialize the InitialMapper object and set up DuckDB connection."""
        self.con = duckdb.connect(db_path)
        self.seq_files = seq_file if isinstance(seq_file, list) else [seq_file]
        self.design_file_path = str(Path(design_file_path).resolve()) if design_file_path else None
        self.bc_objects = bc_objects
        self.reverse_complement = reverse_complement
        self.umi_object = umi_object
        self.step_name = step_name

        self.cols = [bc_object.name for bc_object in bc_objects]
        cols_str = "_".join(self.cols)
        self.table_prefix = f"{step_name}_{cols_str}_"
        self.test_n_reads = None

    @time_it
    def read_fastq(self):
        """
        Read sequences from FASTQ/TXT files and create an initial sequence table.

        If `self.test_n_reads` is set, only a limited number of reads will be read.
        """        
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
                    FROM read_csv_auto('{f_path}', header=False)
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

    def _extract_sequence_object(self, seq_obj):
        """
        Extract a barcode or UMI sequence into a new column and compute quality.

        Args:
            seq_obj (object): Barcode or UMI object with attributes:
                               - name (str)
                               - preceder (str)
                               - post (str)
                               - length (int)
        """
        con = self.con

        # Build regex
        if seq_obj.post == "" or seq_obj.preceder == "":
            regex = f"{seq_obj.preceder}(.{{1,{seq_obj.length}}}){seq_obj.post}"
        else:
            regex = f"{seq_obj.preceder}(.*){seq_obj.post}"
        print(f"Regex for {seq_obj.name}: {regex}")

        # Add sequence column
        con.execute(f"""
            ALTER TABLE seq ADD COLUMN {seq_obj.name} TEXT;
            UPDATE seq
            SET {seq_obj.name} = coalesce(regexp_extract(sequence::VARCHAR, '{regex}', 1), '');
        """)

        # Add quality flag column
        qual_col = f"{seq_obj.name}_qual"
        con.execute(f"ALTER TABLE seq ADD COLUMN {qual_col} BOOLEAN;")
        con.execute(f"""
            UPDATE seq
            SET {qual_col} = LENGTH({seq_obj.name}) = {seq_obj.length};
        """)

    @time_it
    def extract_barcodes(self):
        """Extract all barcode columns and compute their quality flags"""
        print(f"Extracting {len(self.bc_objects)} barcodes...")
        for bc in self.bc_objects:
            self._extract_sequence_object(bc)

    @time_it
    def extract_umi(self):
        """Extract the UMI sequence and compute its quality flag"""
        if self.umi_object:
            print("Extracting UMI...")
            self._extract_sequence_object(self.umi_object)

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

    def _run_pipeline(self):
        """Internal helper: runs the full mapping pipeline."""
        self.read_fastq()
        self.reverse_complement_fastq()
        self.extract_barcodes()
        self.extract_umi()
        self.merge_design()
        print("Mapping complete.")
    
    def create_map(self):
        """Run the full mapping pipeline on all reads."""
        self.test_n_reads = None  # ensure no test limit
        self._run_pipeline()
    
    def create_test_map(self, test_n_reads: int = 100):
        """
        Run the mapping pipeline on a limited number of reads for testing.

        Args:
            test_n_reads (int): Number of reads to process for testing. Defaults to 100.
        """        
        self.test_n_reads = test_n_reads
        self._run_pipeline()
        self.test_n_reads = None  # reset after test   

    def preview_map(self):
        """
        Preview the first 5 rows of the mapped table.

        Returns:
            pandas.DataFrame: A DataFrame of the first 5 rows in the initial mapping table.
        """
        print(f"{self.table_prefix}initial")

        # Count total rows
        total_rows = self.con.execute(f"SELECT COUNT(*) FROM {self.table_prefix}initial").fetchone()[0]
        print(f"Total rows: {total_rows}")
        
        return self.con.execute(f"SELECT * FROM {self.table_prefix}initial LIMIT 5;").df()
