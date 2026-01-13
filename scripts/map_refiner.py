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
import os

from scripts import plotting
from scripts.preprocess import time_it
from scripts.error_correct import run_whitelist_on_concat_domains
from scripts import error_correct

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


    def __init__(self, 
                 db_path, 
                 bc_objects, 
                 column_pairs, 
                 reads_threshold,
                 map_order=None,
                 step_name="", 
                 descriptor="", 
                 design_file = None, 
                 output_figures_path=None, 
                 min_fraction_major_target = 0.9,  
                 manual_ec_threshold = True, 
                 umi_object = None):
        
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.bc_objects = [bc for bc in bc_objects]
        self.cols = [bc.name for bc in self.bc_objects]
        self.column_pairs = column_pairs
        self.reads_threshold = reads_threshold
        self.step_name = step_name
        self.descriptor = descriptor
        self.design_file = design_file
        self.output_figures_path = output_figures_path
        self.min_fraction_major_target = min_fraction_major_target
        self.manual_ec_threshold = manual_ec_threshold
        self.umi_object = umi_object

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

    @time_it
    def create_quality(self, previous_table="initial"):
        """
        Keep only rows where all _qual columns are TRUE.
        """
        print("\nFiltering to high-quality reads...")
    
        # Start with barcode *_qual columns
        qual_columns = [f"{c}_qual = TRUE" for c in self.cols]
        
        # Add UMI *_qual if umi_object exists
        if self.umi_object is not None:
            qual_columns.append(f"{self.umi_object.name}_qual = TRUE")
        
        # Combine into WHERE clause
        where_clause = " AND ".join(qual_columns)
    
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

        if self.output_figures_path:
            sns.set_style('ticks')
            fig, ax = plt.subplots(figsize = (5,2.5), dpi = 300)
            self.plot_reads_histogram(grouped_table_name, ax = ax, edgecolor = 'none', bins = 100)
            plt.title("Grouped " + group_cols_sql)
            
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

        if self.output_figures_path:
            fig, ax = plt.subplots(figsize = (5,2.5), dpi = 300)
            self.plot_reads_histogram(previous_table, ax = ax, edgecolor = 'none', bins = 100)
            ax.axvline(self.reads_threshold, color = 'red')
            
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
        """
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(exist_ok=True)
    
        print(f"\n=== Applying whitelist for {self.step_name} ===")
    
        table_name = prev_table
        new_table = self._prefixed("error_corrected")
    
        # Generate FASTQ for whitelist
        fastq_path = self.generate_fastq_for_whitelist(output_dir, prev_table=table_name)
        fastq_path = Path(fastq_path)
    
        expected_whitelist = output_dir / f"{fastq_path.stem}_whitelist.txt"
        expected_log = output_dir / f"{fastq_path.stem}_whitelist.log"
    
        # Run or reuse whitelist
        if expected_whitelist.exists() and expected_log.exists():
            print(f"Whitelist already exists, skipping:\n  {expected_whitelist}")
            whitelist_path = expected_whitelist
        
        else:
            # Infer expected_bc_count if needed
            if self.manual_ec_threshold:
        
                print("Inferring expected barcode count from grouped reads...")
        
                # Build grouped counts from prev_table
                group_cols_sql = ", ".join(self.cols)
                tmp_grouped = "__tmp_bc_grouped"
        
                self.con.execute(f"""
                    CREATE OR REPLACE TABLE {tmp_grouped} AS
                    SELECT {group_cols_sql}, COUNT(*) AS count
                    FROM {table_name}
                    GROUP BY {group_cols_sql}
                """)
        
                # Prompt for reads threshold if needed (same behavior as create_thresholded)
                if not self.reads_threshold:
                    if self.output_figures_path:
                        # fig, ax = plt.subplots(figsize=(5, 2.5), dpi=300)
                        # self.plot_reads_histogram(tmp_grouped, ax = ax, edgecolor = 'none', bins = 100)
                        fig, ax = plt.subplots(figsize=(5, 2.5), dpi=300)
                        df_tmp = self.con.execute(f"SELECT count FROM {tmp_grouped}").df()
                        sns.histplot(df_tmp["count"], bins=100, ax=ax, log_scale = (True, True))
                        plt.title("Grouped barcode read counts")
                        plt.show()

                        filename = os.path.join(self.output_figures_path, f"{tmp_grouped}.png")
                        plt.savefig(filename, bbox_inches="tight")
        
                    try:
                        user_input = input("Enter minimum read count threshold (default = 5): ")
                        self.reads_threshold = int(user_input) if user_input.strip() != "" else 5
                    except ValueError:
                        print("Invalid input. Using default threshold of 5.")
                        self.reads_threshold = 5
        
                print(f"Using reads threshold of {self.reads_threshold}")
        
                # Count how many groups exceed threshold
                expected_bc_count = self.con.execute(f"""
                    SELECT COUNT(*)
                    FROM {tmp_grouped}
                    WHERE count > {self.reads_threshold}
                """).fetchone()[0]
        
                print(f"Inferred expected barcode count: {expected_bc_count}")
        
                self.con.execute(f"DROP TABLE IF EXISTS {tmp_grouped}")

                print("Running umi_tools whitelist...")
                whitelist_outputs = run_whitelist_on_concat_domains(
                    fastq_path=fastq_path,
                    output_dir=output_dir,
                    prefix=None,
                    set_cell_number=expected_bc_count
                )
                whitelist_path = whitelist_outputs["whitelist"]


            else:
                print(f"Using automatic detection of expected barcode count.")

                print("Running umi_tools whitelist...")
                whitelist_outputs = run_whitelist_on_concat_domains(
                    fastq_path=fastq_path,
                    output_dir=output_dir,
                    prefix=None
                )
                whitelist_path = whitelist_outputs["whitelist"]

        
            
        # Load whitelist mapping
        mapping_df = error_correct.convert_txt_to_whitelist_mapping_df_from_path(whitelist_path)
        print(f"Unique canonical barcodes: {mapping_df['canonical'].nunique()}")
    
        self.con.execute("CREATE OR REPLACE TABLE barcode_map AS SELECT * FROM mapping_df;")
    
        # Join + replace barcodes
        concat_expr = " || ".join(self.cols)
        self.con.execute(f"""
            CREATE OR REPLACE TABLE {new_table} AS
            SELECT t.*, b.canonical AS concat_canonical
            FROM {table_name} AS t
            JOIN barcode_map AS b
            ON {concat_expr} = b.original
        """)
    
        start_idx = 0
        for bc in self.bc_objects:
            self.con.execute(f"""
                UPDATE {new_table}
                SET {bc.name} = SUBSTR(concat_canonical, {start_idx + 1}, {bc.length})
            """)
            self.con.execute(f"""
                UPDATE {new_table}
                SET {bc.name}_qual = LENGTH({bc.name}) = {bc.length}
            """)
            start_idx += bc.length
    
        # Recalculate Designed
        self.con.execute(f"""
            ALTER TABLE {new_table} DROP COLUMN IF EXISTS Designed;
        """)
        self.merge_design(prev_table=new_table)
    
        self.con.execute(f"""
            ALTER TABLE {new_table} DROP COLUMN concat_canonical;
        """)
    
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

    @time_it
    def create_unique_target(self, previous_table):
        """
        For each (key_cols, target_cols) pair, independently keep keys where the most abundant target 
        accounts for ≥ self.min_fraction_major_target of reads. Then, keep rows that pass all filters.
        """
        frac = self.min_fraction_major_target
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

    def refine_map_from_db(self, save_name = None, should_check_exists = False):
        self.should_check_exists = False
        
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

        return plotting.plot_loss_helper(ax=None, 
                                         palette="rocket_r", 
                                         text_offset = 0, 
                                         show_background = True,
                                         default_map_order = self.DEFAULT_MAP_ORDER, 
                                         output_figures_path = self.output_figures_path,
                                         table_prefix_with_descriptor=self.table_prefix_with_descriptor, 
                                         df=df)
        
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
        return plotting.plot_error_correction(self.output_figures_path,
                                              self.table_prefix_with_descriptor,
                                              save_dir,
                                              plot)

    @time_it
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
        return plotting.plot_all_whitelists_from_summary(summary_df, n_cols, dpi)

