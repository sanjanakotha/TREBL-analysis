import duckdb
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white", context="talk")
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
import matplotlib.cm as cm
from pathlib import Path
import plotting
import os
from preprocess import time_it
from error_correct import run_whitelist_on_concat_domains
import error_correct

class MapRefiner:
    """
    A class to refine sequencing maps made through InitialMapper.
    """

    DEFAULT_MAP_ORDER = [
        "initial",
        "quality",
        "error_corrected",   # new step
        "grouped",
        "thresholded",
        "barcode_exists",       
        "unique_target",
        "designed"
    ]


    def __init__(self, db_path, bc_objects, column_pairs, should_check_exists, 
                 design_check, reads_threshold=None, map_order=None, step_name="", 
                 descriptor="", design_file = None, plot_histograms=True, 
                 output_figures_path=None, min_fraction = 0.9, expected_bc_count = None):
        
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.bc_objects = bc_objects
        self.cols = [bc_object.name for bc_object in bc_objects]
        self.column_pairs = column_pairs
        self.design_check = design_check
        self.reads_threshold = reads_threshold
        self.step_name = step_name
        self.descriptor = descriptor
        self.design_file = design_file
        self.plot_histograms = plot_histograms
        self.output_figures_path = output_figures_path
        self.should_check_exists = should_check_exists
        self.min_fraction = min_fraction
        self.expected_bc_count = expected_bc_count

        cols_str = "_".join(self.cols)
        self.table_prefix = f"{step_name}_{cols_str}_"
        
        # --- Build a stable base prefix that never changes across descriptors ---
        # Example: step1_AD_ADBC_RPBC_
        cols_str = "_".join(self.cols)
        if step_name:
            self.table_prefix_base = f"{step_name}_{cols_str}_"
        else:
            self.table_prefix_base = f"{cols_str}_"
    
        # --- Optionally append descriptor later ---
        # This one changes per run if you specify a descriptor
        # e.g. "step1_AD_ADBC_RPBC_strict_"
        if descriptor:
            self.table_prefix_with_descriptor = f"{self.table_prefix_base}{descriptor}_"
        else:
            self.table_prefix_with_descriptor = self.table_prefix_base
    
        print("Base prefix (stable across descriptors):", self.table_prefix_base)
        print("Full prefix for this instance:", self.table_prefix_with_descriptor)
        print()

        valid_steps = [step for step in self.DEFAULT_MAP_ORDER if step != "initial"]

        # Always start with 'initial'
        if map_order is None:
             # --- Map order setup ---
            print("Default map order:")
            for i, name in enumerate(valid_steps, 1):
                print(f"{i}. {name}")
            print()
            
            user_input = input(
                f"Enter map order by numbers for steps after 'initial', comma-separated (e.g., 1,3,2), or press Enter for default: "
            )
            if user_input.strip() == "":
                self.map_order = ["initial"] + valid_steps
            else:
                try:
                    indices = [int(x.strip()) - 1 for x in user_input.split(",")]
                    chosen = [valid_steps[i] for i in indices]
                    self.map_order = ["initial"] + chosen
                except Exception as e:
                    print(f"Invalid input ({e}), using default order.")
                    self.map_order = ["initial"] + valid_steps
        else:
            invalid = [step for step in map_order if step not in valid_steps]
            if invalid:
                raise ValueError(f"Invalid map step(s) in custom order: {invalid}")
            self.map_order = ["initial"] + map_order

        print("Using the following step order:")
        for i, name in enumerate(self.map_order, 1):
            print(f"{i}. {name}")
        print()


    def show_tables(self):
        tables = self.con.execute("SHOW TABLES").fetchall()
        # for table in tables:
        #     print(table)
        return list(tables)
    
    def _prefixed(self, table_name):
        """
        Add prefix to table names.
        - 'initial' always uses base prefix.
        - 'grouped' uses base prefix if previous step in map_order is 'initial',
          otherwise uses descriptor prefix.
        """
        if table_name == "initial":
            return f"{self.table_prefix_base}{table_name}"
    
        elif table_name == "grouped":
            # Determine the index of "grouped" and check what comes before it
            if "grouped" in self.map_order:
                idx = self.map_order.index("grouped")
                prev_step = self.map_order[idx - 1] if idx > 0 else None
            else:
                prev_step = None
    
            # If grouped follows "initial", use base prefix (no descriptor)
            if prev_step == "initial":
                return f"{self.table_prefix_base}initial_grouped"
            else:
                return f"{self.table_prefix_with_descriptor}{table_name}"
    
        else:
            # All other tables use descriptor prefix
            return f"{self.table_prefix_with_descriptor}{table_name}"


    # -------------------------------
    # Table creation functions
    # -------------------------------

    def check_exists(self, check_table_name):
        if self.should_check_exists:
            # Get list of existing tables
            existing_tables = [t[0] for t in self.con.execute("SHOW TABLES").fetchall()]
            
            # Only skip if table exists AND ends with "_initial" or "_initial_grouped"
            if check_table_name in existing_tables and (
                check_table_name.endswith("_initial") or check_table_name.endswith("_initial_grouped") or check_table_name.endswith("_loss_table")
            ):
                print(f"Skipping — table {check_table_name} already exists and is initial/grouped.")
                return True
            else:
                return False
        else:
            return False
            
    @time_it
    def create_initial(self, parquet_path):
        print(f"Reading initial map from {parquet_path} as {self._prefixed('initial')}...")
        if not self.check_exists(self._prefixed('initial')):
            self.con.execute(f"""
                CREATE OR REPLACE TABLE {self._prefixed('initial')} AS
                SELECT *
                FROM read_parquet('{parquet_path}')
            """)

    # @time_it
    # def create_quality_designed(self, previous_table="initial"):
    #     print("Filtering to quality and designed...")
        
    #     # Keep only rows where all *_qual columns are TRUE
    #     where_clause = " AND ".join([f"{c}_qual = TRUE" for c in self.cols])
    
    #     # Require Designed == 1 (string or int safe)
    #     designed_filter = " AND CAST(Designed AS INTEGER) = 1" if self.design_check else ""
    
    #     self.con.execute(f"""
    #         CREATE OR REPLACE TABLE {self._prefixed('quality_designed')} AS
    #         SELECT *
    #         FROM {self._prefixed(previous_table)}
    #         WHERE {where_clause}
    #         {designed_filter}
    #     """)

    @time_it
    def create_quality(self, previous_table="initial"):
        """
        Keep only rows where all *_qual columns are TRUE.
        """
        print("\nFiltering to high-quality reads...")
    
        # Build WHERE clause using the known barcode names
        where_clause = " AND ".join([f"{c}_qual = TRUE" for c in self.cols])
    
        self.con.execute(f"""
            CREATE OR REPLACE TABLE {self._prefixed('quality')} AS
            SELECT *
            FROM {self._prefixed(previous_table)}
            WHERE {where_clause};
        """)
    
        print(f"Created table: {self._prefixed('quality')} — filtered for TRUE in all *_qual columns.")
    
    @time_it
    def create_designed(self, previous_table="quality"):
        """
        Keep only rows where Designed == 1.
        """
        print("\nFiltering to designed sequences...")
    
        self.con.execute(f"""
            CREATE OR REPLACE TABLE {self._prefixed('designed')} AS
            SELECT *
            FROM {self._prefixed(previous_table)}
            WHERE CAST(Designed AS INTEGER) = 1;
        """)
    
        print(f"Created table: {self._prefixed('designed')} — kept only Designed == 1.")

    @time_it
    def barcode_exists(self, previous_table="initial"):
        """
        Remove rows where any barcode column (except AD) is NULL or an empty string.
        """
        print("Removing rows with null or empty barcodes (excluding AD)...")
        prev_table = self._prefixed(previous_table)
        output_table = self._prefixed("barcode_exists")
    
        # Filter columns: exclude 'AD'
        barcode_cols = [c for c in self.cols if c != "AD"]
    
        if barcode_cols:
            # Build condition to detect NULL or empty strings in barcode columns
            null_empty_condition = " OR ".join([f"{c} IS NULL OR TRIM({c}) = ''" for c in barcode_cols])
            where_clause = f"WHERE NOT ({null_empty_condition})"
        else:
            # No barcode columns to filter
            where_clause = ""
    
        self.con.execute(f"""
            CREATE OR REPLACE TABLE {output_table} AS
            SELECT *
            FROM {prev_table}
            {where_clause}
        """)
    
    @time_it
    def create_grouped(self, previous_table="quality_designed"):
        print(f"Grouping {self._prefixed(previous_table)}...")

        grouped_table_name = self._prefixed('grouped')
        group_cols_sql = ", ".join(self.cols)
        
        if not self.check_exists(grouped_table_name):
            qual_cols = [f"{c}_qual" for c in self.cols]
            qual_cols_sql = ", ".join(qual_cols + ["Designed"])
            self.con.execute(f"""
                CREATE OR REPLACE TABLE {grouped_table_name} AS
                SELECT {group_cols_sql},
                    COUNT(*) AS count,
                    {qual_cols_sql}
                FROM {self._prefixed(previous_table)}
                GROUP BY {group_cols_sql}, {qual_cols_sql}
                ORDER BY count DESC
            """)

        if self.plot_histograms:
            sns.set_style('ticks')
            fig, ax = plt.subplots(figsize = (5,2.5), dpi = 300)
            self.plot_reads_histogram(grouped_table_name, ax = ax, edgecolor = 'none', bins = 100)
            plt.title("Grouped " + group_cols_sql)
            
            if self.output_figures_path:
                grouped_table_name = self._prefixed('grouped')
                filename = os.path.join(self.output_figures_path, f"{grouped_table_name}.png")
                plt.savefig(filename, bbox_inches="tight")

            plt.show()

    @time_it
    def create_thresholded(self, previous_table):
        print("Thresholding...")
        
        if not self.reads_threshold:    
            try:
                user_input = input("Enter minimum read count threshold (default = 5): ")
                self.reads_threshold = int(user_input) if user_input.strip() != "" else 5
            except ValueError:
                print("Invalid input. Using default threshold of 5.")
                self.reads_threshold = 5
        
        print("Using reads threshold of " + str(self.reads_threshold) + ".")

        if self.plot_histograms:
            fig, ax = plt.subplots(figsize = (5,2.5), dpi = 300)
            self.plot_reads_histogram(previous_table, ax = ax, edgecolor = 'none', bins = 100)
            ax.axvline(self.reads_threshold, color = 'red')

            if self.output_figures_path:
                grouped_table_name = self._prefixed('grouped')
                filename = os.path.join(self.output_figures_path, f"{self._prefixed('thresholded')}.png")
                plt.savefig(filename, bbox_inches="tight")

            plt.show()

        self.con.execute(f"""
            CREATE OR REPLACE TABLE {self._prefixed('thresholded')} AS
            SELECT *
            FROM {self._prefixed(previous_table)}
            WHERE count > {self.reads_threshold}
        """)

    # ERR CORRECT
    @time_it
    def generate_fastq_for_whitelist(self, output_dir=None, suffix="_barcodes_extracted.fastq", prev_table = "quality"):
        """
        Generate a synthetic FASTQ file from the filtered previous table, using
        concatenated barcode for reads that pass quality checks.
    
        Args:
            output_dir (str or Path, optional): Directory to save the FASTQ. Defaults to current working directory.
            suffix (str): Suffix for FASTQ filename.
            prev_table (str, optional): Name of the previous DuckDB table to use instead of initial table.
        
        Returns:
            str: Path to the generated FASTQ file.
        """
        con = self.con
    
        # Resolve output directory
        if output_dir is None:
            output_dir = Path(".")
        else:
            output_dir = Path(output_dir).resolve()
            output_dir.mkdir(exist_ok=True)
    
        # Construct output FASTQ path
        output_fastq = output_dir / f"{self.table_prefix}filtered{suffix}"
        print(f"Generating FASTQ: {output_fastq}")
    
        # Get all barcode and quality column names
        bc_cols = self.cols
        qual_cols = [f"{name}_qual" for name in bc_cols]
    
        # Build quality filter (require all barcodes to pass length check)
        qual_filter = " AND ".join([f"{qc} = TRUE" for qc in qual_cols])
    
        # Construct concatenation expression for barcodes (+ optional UMI)
        concat_expr = " || ".join(bc_cols)
    
        # Use previous table
        table_name = prev_table 
    
        # Query the concatenated reads that pass quality checks
        query = f"""
            SELECT {concat_expr} AS sequence
            FROM {table_name}
            WHERE {qual_filter};
        """
        seq_df = con.execute(query).df()
    
        # Construct FASTQ records
        with open(output_fastq, "w") as f:
            for i, seq in enumerate(seq_df["sequence"], start=1):
                # Dummy FASTQ format: 4 lines per record
                # Adding "A" as dummy UMI
                f.write(f"@read_{i}\nA{seq}\n+\n{'I' * (len(seq) + 1)}\n")
    
        print(f"Wrote {len(seq_df)} reads to {output_fastq}")
        return str(output_fastq)
    
    
    def apply_whitelist(self, output_dir, prev_table=None):
        """
        Apply umi_tools whitelist to previous table: replace barcodes and AD sequences with canonical ones,
        then recalculate quality flags and 'Designed' based on canonical barcode.
        Saves results to a new table prefixed with error_corrected.
    
        Args:
            output_dir (str or Path): Directory where FASTQ / outputs will be saved.
            prev_table (str, optional): Name of previous DuckDB table to use instead of initial.
        """
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(exist_ok=True)
    
        print(f"\n=== Applying whitelist for {self.step_name} ===")
    
        # Use previous table or default initial table
        table_name = prev_table #or f"{self.table_prefix}initial"
        new_table = self._prefixed("error_corrected")
    
        # Generate FASTQ for whitelist
        fastq_path = self.generate_fastq_for_whitelist(output_dir, prev_table=table_name)
        fastq_path = Path(fastq_path)
    
        # Run umi_tools whitelist
        expected_whitelist = output_dir / f"{fastq_path.stem}_whitelist.txt"
        expected_log = output_dir / f"{fastq_path.stem}_whitelist.log"
               
        if expected_whitelist.exists() and expected_log.exists():
            print(f"Whitelist already exists, skipping:\n  {expected_whitelist}")
            whitelist_path = expected_whitelist
        elif self.expected_bc_count:
            print("Running umi_tools whitelist...")
            whitelist_outputs = run_whitelist_on_concat_domains(
                fastq_path=fastq_path,
                output_dir=output_dir,
                prefix=None,
                set_cell_number = self.expected_bc_count
            )
        else:
            print("Running umi_tools whitelist...")
            whitelist_outputs = run_whitelist_on_concat_domains(
                fastq_path=fastq_path,
                output_dir=output_dir,
                prefix=None
            )
            whitelist_path = whitelist_outputs["whitelist"]
    
       # Step 3: load mapping_df
        mapping_df = error_correct.convert_txt_to_whitelist_mapping_df_from_path(whitelist_path)
        print(f"Unique canonical barcodes: {len(mapping_df['canonical'].unique())}")
    
        # Step 4: push mapping_df to DuckDB
        self.con.execute("CREATE OR REPLACE TABLE barcode_map AS SELECT * FROM mapping_df;")
    
        # Step 5: join prev_table to mapping_df on concatenated barcodes (filters unmatched reads)
        concat_expr = " || ".join(self.cols)
        self.con.execute(f"""
            CREATE OR REPLACE TABLE {new_table} AS
            SELECT t.*, b.canonical AS concat_canonical
            FROM {table_name} AS t
            JOIN barcode_map AS b
            ON {concat_expr} = b.original
        """)
    
        # Step 6: update individual barcode columns with canonical subsequences
        start_idx = 0
        for bc in self.bc_objects:
            end_idx = start_idx + bc.length
            self.con.execute(f"""
                UPDATE {new_table}
                SET {bc.name} = SUBSTR(concat_canonical, {start_idx + 1}, {bc.length})
            """)
            # recalc quality
            self.con.execute(f"""
                UPDATE {new_table}
                SET {bc.name}_qual = LENGTH({bc.name}) = {bc.length}
            """)
            start_idx = end_idx

        # Step 7a: drop existing Designed column if it exists
        self.con.execute(f"""
            ALTER TABLE {new_table} DROP COLUMN IF EXISTS Designed;
        """)
    

        # Step 7: recalc Designed using canonical AD column
        self.merge_design(prev_table=new_table)
    
        # Step 8: drop temporary concat_canonical
        self.con.execute(f"ALTER TABLE {new_table} DROP COLUMN concat_canonical;")
    
        print(f"Whitelist application complete for {self.step_name} at {new_table}")

    @time_it
    def merge_design(self, prev_table=None):
        """
        Merge with design file if provided, or create default Designed column.
    
        Args:
            prev_table (str, optional): Name of the DuckDB table to merge design info into.
                                         Defaults to `seq` if not provided.
        """
        con = self.con
        table_to_use = prev_table or "seq"
        target_table = table_to_use  # output table is now the same as input table
    
        if self.design_file:
            print("Merging with design file...")
            # Load design file
            con.execute(f"""
                CREATE OR REPLACE TABLE design AS
                SELECT CAST(column0 AS VARCHAR) AS AD
                FROM read_csv_auto('{self.design_file}', header=False)
            """)
            # Merge design info into the table
            con.execute(f"""
                CREATE OR REPLACE TABLE {target_table} AS
                SELECT t.*, CASE WHEN d.AD IS NOT NULL THEN 1 ELSE 0 END AS Designed
                FROM {table_to_use} AS t
                LEFT JOIN design AS d USING(AD);
            """)
        else:
            # No design file, default Designed = 1
            con.execute(f"""
                CREATE OR REPLACE TABLE {target_table} AS
                SELECT *, 1 AS Designed
                FROM {table_to_use};
            """)
    
        # Drop sequence column if present
        columns = [col for col in con.execute(f"PRAGMA table_info('{target_table}')").fetchall()]
        if any(c[1] == 'sequence' for c in columns):
            con.execute(f"ALTER TABLE {target_table} DROP COLUMN sequence;")
    
    @time_it
    def create_error_corrected(self, previous_table="quality"):
        """
        Perform error correction by generating a FASTQ from grouped barcodes,
        applying umi_tools whitelist, and updating the database with corrected barcodes.
        """
        print(f"\n=== Running error correction step on {self._prefixed(previous_table)} ===")
        output_dir = self.output_figures_path or "./error_corrected"
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate FASTQ + run whitelist-based correction
        try:
            self.apply_whitelist(output_dir, prev_table = self._prefixed(previous_table))
            
        except Exception as e:
            print(f"⚠️ Error correction failed: {e}")


    # @time_it
    # def create_unique_target(self, previous_table):
    #     print("Filtering keys to keep only those mapping to a single target...")
    
    #     previous_table = self._prefixed(previous_table)
    #     tmp_tables = []
    
    #     for i, (key_cols, target_cols) in enumerate(self.column_pairs):
    #         tmp_table = self._prefixed(f"tmp_map4_{i}")
    #         tmp_tables.append(tmp_table)
    
    #         # Single or multi-column key
    #         if isinstance(key_cols, str):
    #             key_expr = f"{key_cols}"
    #         else:
    #             key_expr = " || '-' || ".join([f"{c}" for c in key_cols])
    
    #         # Single or multi-column target
    #         if isinstance(target_cols, str):
    #             target_expr = f"{target_cols}"
    #         else:
    #             target_expr = " || '-' || ".join([f"{c}" for c in target_cols])

    #         print(f"\tChecking each {key_expr} only maps to one {target_expr}.")
            
    #         query = f"""
    #              CREATE OR REPLACE TABLE {tmp_table} AS
    #             SELECT *
    #             FROM {previous_table}
    #             WHERE {key_expr} IN (
    #                 SELECT {key_expr}
    #                 FROM {previous_table}
    #                 GROUP BY {key_expr}
    #                 HAVING COUNT(DISTINCT {target_expr}) = 1
    #             )
    #         """            
    #         # Only keep keys that map to exactly one target
    #         self.con.execute(query)
    
    #         previous_table = tmp_table
    
    #     # Finalize
    #     self.con.execute(f"""
    #         CREATE OR REPLACE TABLE {self._prefixed('unique_target')} AS
    #         SELECT * FROM {previous_table};
    #     """)
    
    #     # Cleanup
    #     for t in tmp_tables:
    #         self.con.execute(f"DROP TABLE IF EXISTS {t};")
    

    @time_it
    def create_unique_target(self, previous_table):
        """
        For each (key_cols, target_cols) pair, independently keep keys where the most abundant target 
        accounts for ≥ self.min_fraction of reads. Then, keep rows that pass all filters.
        """
        frac = self.min_fraction
        print(f"Filtering so that at least {frac*100:.0f}% of reads come from the most abundant target...")
    
        base_table = self._prefixed(previous_table)
        filter_tables = []
    
        # Step 1: build independent filters for each (key, target) mapping
        for i, (key_cols, target_cols) in enumerate(self.column_pairs):
            tmp_name = self._prefixed(f"unique_target_step_{i}")
            filter_tables.append(tmp_name)
    
            # Build key/target expressions
            key_expr = (
                key_cols if isinstance(key_cols, str)
                else " || '-' || ".join([f"{c}" for c in key_cols])
            )
            target_expr = (
                target_cols if isinstance(target_cols, str)
                else " || '-' || ".join([f"{c}" for c in target_cols])
            )
    
            print(f"\tProcessing mapping {i+1}: {key_expr} → {target_expr}")
    
            query = f"""
                CREATE OR REPLACE TABLE {tmp_name} AS
                WITH counts AS (
                    SELECT 
                        {key_expr} AS computed_key,
                        {target_expr} AS computed_target,
                        SUM(count) AS cnt
                    FROM {base_table}
                    GROUP BY computed_key, computed_target
                ),
                totals AS (
                    SELECT 
                        computed_key,
                        SUM(cnt) AS total_cnt,
                        MAX(cnt) AS max_cnt
                    FROM counts
                    GROUP BY computed_key
                )
                SELECT computed_key
                FROM totals
                WHERE CAST(max_cnt AS FLOAT) / total_cnt >= {frac};
            """
    
            self.con.execute(query)
    
        # Step 2: intersect all filter tables
        final_name = self._prefixed("unique_target")
        final_conditions = []
    
        for tbl, (key_cols, _) in zip(filter_tables, self.column_pairs):
            key_expr = (
                key_cols if isinstance(key_cols, str)
                else " || '-' || ".join([f"{c}" for c in key_cols])
            )
            final_conditions.append(f"{key_expr} IN (SELECT computed_key FROM {tbl})")
    
        query_final = f"""
            CREATE OR REPLACE TABLE {final_name} AS
            SELECT *
            FROM {base_table} AS base
            WHERE {" AND ".join(final_conditions)};
        """
    
        self.con.execute(query_final)
    
        # Step 3: optional cleanup of intermediate tables
        for tbl in filter_tables:
            self.con.execute(f"DROP TABLE IF EXISTS {tbl};")
    
        print(f"Created filtered table: {final_name}")
        return final_name


    # -------------------------------
    # Pipelines
    # -------------------------------
    @time_it
    def load_csv(self, csv_path, table_name):
        """
        Load a CSV file as the 'initial' table in the database.
    
        Parameters
        ----------
        csv_path : str
            Path to the CSV file.
        """
        table_name = self._prefixed(table_name)
        print(f"Reading CSV from {csv_path} into {table_name}...")
    
        if not self.check_exists(table_name):
            self.con.execute(f"""
                CREATE OR REPLACE TABLE {table_name} AS
                SELECT *
                FROM read_csv_auto('{csv_path}')
            """)
        
        
    # def refine_map_from_parquet(self, parquet_path, save_name = None):
    #     table_to_func = {
    #         "initial": self.create_initial,
    #         "quality_designed": self.create_quality_designed,
    #         "barcode_exists": self.barcode_exists,    # new step
    #         "grouped": self.create_grouped,
    #         "thresholded": self.create_thresholded,
    #         "unique_target": self.create_unique_target
    #     }
    #     prev_table = None
    #     for table_name in self.map_order:
    #         func = table_to_func.get(table_name)
    #         if func is None:
    #             continue
    #         if table_name == "initial":
    #             func(parquet_path)
    #             prev_table = "initial"
    #         else:
    #             func(previous_table=prev_table)
    #             prev_table = table_name

    #     if save_name:
    #         self.con.execute(f"CREATE OR REPLACE TABLE {save_name} AS SELECT * FROM {self._prefixed(prev_table)}")
    #     print("Done.")

    def refine_map_from_db(self, save_name = None):
        table_to_func = {
            "initial": None,
            "quality": self.create_quality,
            "designed": self.create_designed,
            "barcode_exists": self.barcode_exists,
            "grouped": self.create_grouped,
            "error_corrected": self.create_error_corrected,  # new
            "thresholded": self.create_thresholded,
            "unique_target": self.create_unique_target
        }

        prev_table = "initial"
        for table_name in self.map_order:
            func = table_to_func.get(table_name)
            if func is not None:
                func(previous_table=prev_table)
            prev_table = table_name

        if save_name:
            self.con.execute(f"CREATE OR REPLACE TABLE {save_name} AS SELECT * FROM {self._prefixed(prev_table)}")
        print("Done.")
            
    # -------------------------------
    # Table access / saving
    # -------------------------------
    def get_map_df(self, map_name, insert_prefix = True):
        existing_tables = [t[0] for t in self.con.execute("SHOW TABLES").fetchall()]
        prefixed_name = self._prefixed(map_name)
        if prefixed_name not in existing_tables:
            return self.con.execute(f"SELECT * FROM {map_name}").df()
        else:
            return self.con.execute(f"SELECT * FROM {prefixed_name}").df()

    def save_map(self, map_name, map_path):
        self.con.execute(f"COPY {self._prefixed(map_name)} TO '{map_path}' WITH (HEADER, DELIMITER ',')")

    @time_it
    def save_loss_table(self, output_csv_path=None):
        """
        Generates a table summarizing the counts and losses at each step of the mapping process.
    
        Parameters
        ----------
        output_csv_path : str or None
            If provided, saves the resulting DataFrame as a CSV file at this path.
    
        Returns
        -------
        pd.DataFrame
            DataFrame summarizing unique counts, unique AD counts, total reads,
            and percentage losses at each mapping step.
        """
        # Mapping step descriptions
        map_info_dict = {
            "initial": "Initial combinations",
            "grouped": "Grouped counts",
            "thresholded": f"Filtered by # reads > {self.reads_threshold}",
            "unique_target": "Filtered for unique targets",
            "quality": "After quality filtering",
            "designed": "After designed filtering",

        }
    
        lengths = []
        existing_tables = [t[0] for t in self.con.execute("SHOW TABLES").fetchall()]
    
        for table_name in self.map_order:
            prefixed_name = self._prefixed(table_name)
            if prefixed_name in existing_tables:
                # Get unique count
                unique_count = self.con.execute(f"SELECT COUNT(*) FROM {prefixed_name}").fetchone()[0]
                
                # Get column info
                columns = [col[0] for col in self.con.execute(f"DESCRIBE {prefixed_name}").fetchall()]
                
                # Total reads
                total_reads = self.con.execute(f"SELECT SUM(count) FROM {prefixed_name}").fetchone()[0] if "count" in columns else unique_count
                
                # Unique AD count
                unique_AD_count = self.con.execute(f"SELECT COUNT(DISTINCT AD) FROM {prefixed_name}").fetchone()[0] if "AD" in columns else 0
            else:
                # Table missing
                unique_count = None
                total_reads = None
                unique_AD_count = None
    
            # Append info for this step
            description = map_info_dict.get(table_name, table_name)
            lengths.append({
                "map": table_name,
                "description": description,
                "unique_count": unique_count,
                "unique_AD_count": unique_AD_count,
                "total_reads": total_reads
            })
    
        # Build DataFrame
        df = pd.DataFrame(lengths)
    
        # Compute percentages
        prev_count = None
        map1_count = df.loc[df['map'] == 'initial', 'unique_count'].values[0] if not df.empty else None
        percent_prev, percent_map1 = [], []
    
        for count in df['unique_count']:
            # % vs previous step
            percent_prev.append(round(100 * count / prev_count, 2) if count is not None and prev_count is not None else None)
            # % vs initial
            percent_map1.append(round(100 * count / map1_count, 2) if count is not None and map1_count is not None else None)
            prev_count = count if count is not None else prev_count
    
        df['% of previous step'] = pd.Series(percent_prev)
        df['% of total reads'] = pd.Series(percent_map1)
    
        # Fill NaNs in percentage columns only
        df['% of previous step'] = df['% of previous step'].fillna(100)
        df['% of total reads'] = df['% of total reads'].fillna(100)
    
        # Save CSV if requested
        if output_csv_path:
            df.to_csv(output_csv_path, index=False)

        
        # Save summary table to SQL (overwrite if exists)
        loss_table_name = self._prefixed("loss_summary")
        self.con.execute(f"DROP TABLE IF EXISTS {loss_table_name}")
        self.con.execute(f"CREATE TABLE {loss_table_name} AS SELECT * FROM df")
    
        print(f"Saved loss summary table as '{loss_table_name}'")
        
        return df

    # -------------------------------
    # Plotting functions
    # -------------------------------
    def plot_loss(self, ax=None, palette="rocket_r", text_offset = 0, show_background = True):
        """
        Plot loss across mapping steps showing total reads vs unique counts.
        Adds numeric count labels at the end of each bar.
    
        Args:
            save_path (str or None): Optional path to save the plot.
            ax (matplotlib.axes.Axes or None): Optional existing axis to draw on.
            palette (str): Seaborn color palette name.
    
        Returns:
            matplotlib.axes.Axes: The plotted axis.
        """
        # Name of saved summary table
        loss_table_name = self._prefixed("loss_summary")
    
        # Use check_exists() and get_map_df() to get the DataFrame
        if self.check_exists(loss_table_name):
            df = self.get_map_df("loss_summary")
        else:
            df = self.save_loss_table()

        #df = df[df["map"] != "initial"]
    
        # Initialize plot
        if ax is None:
            sns.set(style="white", context="talk")
            fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    
        # Create color mapping
        cmap = sns.color_palette(palette, n_colors=len(self.DEFAULT_MAP_ORDER))
        palette_dict = {name: cmap[i] for i, name in enumerate(self.DEFAULT_MAP_ORDER)}
        colors = [palette_dict.get(name, (0.5, 0.5, 0.5)) for name in df["map"]]
    
        # Bar plots for total reads (light) and unique counts (solid)
        if show_background:
            background_alpha = 0.5
        else:
            background_alpha = 0
        sns.barplot(x="total_reads", y="description", data=df, ax=ax, palette=colors, alpha=background_alpha)
        sns.barplot(x="unique_count", y="description", data=df, ax=ax, palette=colors, alpha = 1)
        #sns.scatterplot(x="unique_AD_count", y="description", data=df, ax=ax, palette=colors, zorder = 20)

        # Draw black line marking unique_AD_count — only show top edge
        # sns.barplot(
        #     x="unique_AD_count", y="description", data=df,
        #     ax=ax, color='none', zorder=20
        # )
            
        # Add count labels to bars
        for i, row in df.iterrows():
            # Unique counts (front bar)
            if pd.notna(row["unique_count"]):
                if row["unique_count"] != row["total_reads"]:
                    ax.text(
                        row["unique_count"] * 1.01,
                        i-text_offset,
                        f"{int(row['unique_count']):,}",
                        va="center",
                        ha="left",
                        fontsize='small',
                        color="black"#, zorder = 20
                    )
            if show_background:
                background_color = 'black'
            else:
                background_color = 'white'
            # Total reads (back bar)
            if pd.notna(row["total_reads"]):
                ax.text(
                    row["total_reads"] * 1.01,
                    i+text_offset,
                    f"{int(row['total_reads']):,}",
                    va="center",
                    ha="left",
                    fontsize='small',
                    color=background_color#, zorder = 20
                )

        # if self.design_file:
        #     design_file_df = pd.read_csv(self.design_file)
        #     ax.axvline(len(design_file_df), color = 'gray', alpha = 0.5, zorder = 0, linestyle = 'dashed')
        #     print(len(design_file_df))
        
        ax.set_xlabel("Read Count (Unique, Total)")
        ax.set_ylabel("Map Step")
        ax.set_yticklabels(df["map"])
        sns.despine(bottom=True)
        ax.set_xticks([])
    
        if self.output_figures_path:
            grouped_table_name = self._prefixed('grouped')
            filename = os.path.join(self.output_figures_path, f"{self.table_prefix_with_descriptor}loss.png")
            plt.savefig(filename, bbox_inches="tight")

        return ax
        
    def plot_reads_histogram(self, previous_table, save_path=None, ax=None, **kwargs):
        map_df = self.get_map_df(previous_table)
        return plotting.plot_reads_histogram(map_df, save_path=save_path, ax=ax, **kwargs)

    def close_connection(self):
        if self.con:
            self.con.close()
            self.con = None

    @time_it
    def plot_error_correction(self, save_dir=None, plot=True):
        """
        Summarize counts of canonical subsequences (barcodes) per AD sequence
        from whitelist files, produce a summary table and optional plots.
    
        Args:
            save_dir (str or Path, optional): Directory to save outputs.
            plot (bool): Whether to produce plots.
    
        Returns:
            pd.DataFrame: summary table of barcodes with counts and group sizes.
        """
        import glob
    
        whitelist_dir = Path(self.output_figures_path or "./error_corrected")
        whitelist_files = sorted(glob.glob(str(whitelist_dir / "*_whitelist.txt")))
        print(whitelist_files)
        
        if not whitelist_files:
            raise FileNotFoundError(f"No whitelist files found in {whitelist_dir}")
    
    
        wl_df = pd.read_csv(whitelist_files[0], sep="\t", header=None,
                            names=["canonical", "collapsed", "largest_count", "counts"])

        # Process collapsed sequences
        wl_df["collapsed_list"] = wl_df["collapsed"].fillna("").apply(lambda x: x.split(",") if x else [])
        wl_df["num_merged"] = wl_df["collapsed_list"].apply(len)

        # Compute rest counts
        def parse_counts(row):
            if pd.isna(row["counts"]) or str(row["counts"]).strip() == "":
                return row["largest_count"], 0
            rest = sum(int(c) for c in str(row["counts"]).split(",") if c.strip())
            largest = row["largest_count"] - rest
            return largest, rest

        wl_df[["largest_count", "rest_count"]] = wl_df.apply(parse_counts, axis=1, result_type="expand")

        summary_df = wl_df
    
        # # Save summary
        # if save_dir:
        #     save_dir = Path(save_dir)
        #     save_dir.mkdir(exist_ok=True)
        #     csv_path = save_dir / f"{self.table_prefix_with_descriptor}whitelist_summary.csv"
        #     summary_df.to_csv(csv_path, index=False)
        #     print(f"Saved whitelist summary: {csv_path}")
    
        #Plot distributions
        if plot:
            sns.set(style="white", context="talk")

            fig = self.plot_all_whitelists_from_summary(summary_df)
            if self.output_figures_path:
                plot_path = Path(self.output_figures_path) / f"{self.table_prefix_with_descriptor}whitelist_summary.png"
                fig.savefig(plot_path, bbox_inches="tight")
                print(f"Saved barcode whitelist plots: {plot_path}")
            plt.show()
    
    def plot_all_whitelists_from_summary(self, summary_df, n_cols=4, dpi=300):
        """
        Plot summaries of multiple barcodes using precomputed summary_df
        instead of reading whitelist files again.
    
        Parameters
        ----------
        summary_df : pd.DataFrame
            DataFrame returned from summarize_canonical_subsequences, must contain
            columns ['barcode', 'canonical', 'num_merged', 'largest_count', 'rest_count'].
        n_cols : int, optional
            Number of panels per row. Default is 4.
        dpi : int, optional
            Figure DPI. Default is 300.
    
        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure containing all barcode summaries.
        """
        n_barcodes = 1
        n_rows = n_barcodes
        n_panels = 4  # same as summarize_whitelist
        
        fig, axs = plt.subplots(n_rows, n_panels, figsize=(5 * n_panels, 5 * n_rows), dpi=dpi)
    
        # Ensure axs is 2D array for consistent indexing
        if n_rows == 1:
            axs = axs.reshape(1, -1)

        i = 0
        df_bc = summary_df#[summary_df["barcode"] == bc_name]
        # df_bc = df_bc[["barcode", "num_merged", "rest_count", "largest_count"]]
        # df_bc = df_bc.groupby("barcode").sum()
        # df_bc.columns = ["barcode", "num_merged", "rest_count", "largest_count"]

        # Panel 1: Bar chart of original vs canonical sequences
        total_original_seqs = df_bc["num_merged"].sum() + len(df_bc)
        total_canonical_seqs = len(df_bc)
        axs[i, 0].bar(["Before", "After"], [total_original_seqs, total_canonical_seqs],
                      color=[sns.color_palette('Paired')[0], sns.color_palette('Paired')[1]])
        axs[i, 0].set_ylabel("Number of sequences")
        for x, y in zip(["Before", "After"], [total_original_seqs, total_canonical_seqs]):
            axs[i, 0].text(x, y + max([total_original_seqs, total_canonical_seqs])*0.02, f"{y:,}", 
                           ha='center', va='bottom', fontsize='medium', weight='bold')

        # Panel 2: Histogram of group sizes
        axs[i, 1].hist(df_bc["num_merged"] + 1, bins=30, edgecolor='black')
        axs[i, 1].set_xlabel("Group size")
        axs[i, 1].set_ylabel("Frequency")

        # Panel 3: Scatter - largest member vs group size
        axs[i, 2].scatter(df_bc["largest_count"], df_bc["num_merged"] + 1, alpha=0.6)
        axs[i, 2].set_xlabel("Reads of largest")
        axs[i, 2].set_ylabel("Group size")

        # Panel 4: Scatter - largest member vs sum of smaller members
        axs[i, 3].scatter(df_bc["largest_count"], df_bc["rest_count"], alpha=0.6)
        axs[i, 3].set_xlabel("Reads of largest")
        axs[i, 3].set_ylabel("Summed merged reads")
        axs[i, 3].set_xscale('log')
        axs[i, 3].set_yscale('log')

        # Plot y=x line for panel 4
        xlims = axs[i, 3].get_xlim()
        ylims = axs[i, 3].get_ylim()
        line_min = min(xlims[0], ylims[0])
        line_max = max(xlims[1], ylims[1])
        axs[i, 3].plot([line_min, line_max], [line_min, line_max], 'r--', alpha=0.7)
        axs[i, 3].set_xlim(xlims)

        # Add row-level label for the barcode
        fig.text(
            -0.02,
            (n_rows - i - 0.5) / n_rows,
            f"Concat",
            fontsize='large',
            weight='bold',
            rotation=90,
            va='center'
        )

        sns.despine()
        plt.tight_layout(pad=2)
        return fig

