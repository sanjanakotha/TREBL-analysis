import duckdb
import pandas as pd

class MapRefiner:
    """
    A class to refine sequencing maps made through BarcodeMapper.

    Args:
        db_path (str): Path for new or existing DuckDB database.
        cols (list[str]): List of barcode or tile column names to process. At least one must be "AD".
        reads_threshold (int): Minimum reads threshold for filtering map4.
        column_pairs (list[tuple]): List of (key_cols, target_cols) for enforcing each target maps to one key.
        key_cols (str or tuple): Column(s) forming the key.
        target_cols (str or tuple): Column(s) forming the target.

    Example:
        >>> refiner = MapRefiner(
        ...     db_path="my_database.duckdb",
        ...     cols=["AD", "AD_BC", "RP_BC"],
        ...     reads_threshold=5,
        ...     column_pairs=[("AD", "RP_BC"), (("HawkBCs", "AD_BC"), "RP_BC")]
        ... )
        >>> refiner.create_map1("reads.parquet")
        >>> refiner.create_map1("folder/*.parquet")
    """

    def __init__(self, db_path, cols, reads_threshold, column_pairs):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.cols = cols
        self.reads_threshold = reads_threshold
        self.column_pairs = column_pairs

    def create_map1(self, parquet_path):
        """
        Load raw sequencing data from Parquet file or pattern into map1.

        Args:
            parquet_path (str): Path to input Parquet file(s).

        Example:
            >>> refiner.create_map1("reads.parquet")
            >>> refiner.create_map1("folder/*.parquet")
        """
        self.con.execute(f"""
            CREATE OR REPLACE TABLE map1 AS
            SELECT *
            FROM read_parquet('{parquet_path}')
        """)

    def create_map2(self):
        """
        Remove low-quality or non-designed rows from map1 to create map2.

        Specifically:
            1. Removes low-quality barcodes or tiles.
            2. Removes non-designed tiles.

        Example:
            >>> refiner.create_map2()
        """
        cols_to_keep = ", ".join(self.cols)
        where_clause = " AND ".join([f"{c}_qual NOT IN (0, FALSE)" for c in self.cols])

        self.con.execute(f"""
            CREATE OR REPLACE TABLE map2 AS
            SELECT {cols_to_keep}
            FROM map1
            WHERE {where_clause}
            AND Designed NOT IN (0, FALSE)
        """)

    def create_map3(self):
        """
        Group map2 by barcode/tile columns and count occurrences to create map3.

        Example:
            >>> refiner.create_map3()
        """
        group_cols_sql = ", ".join(self.cols)

        self.con.execute(f"""
            CREATE OR REPLACE TABLE map3 AS
            SELECT {group_cols_sql}, COUNT(*) AS count
            FROM map2
            GROUP BY {group_cols_sql}
            ORDER BY count DESC
        """)

    def create_map4(self):
        """
        Filter map3 by read count threshold to create map4.

        Example:
            >>> refiner.create_map4()
        """
        self.con.execute(f"""
            CREATE OR REPLACE TABLE map4 AS
            SELECT *
            FROM map3
            WHERE count > {self.reads_threshold}
        """)

    def create_map5(self):
        """
        Enforce mapping constraints from column_pairs on map4 to create map5.

        Each target combination should map to exactly one key combination.
        Supports single-column or multi-column keys and targets.

        Example:
            >>> refiner.create_map5()
        """
        current_table = "map4"

        for i, (key_cols, target_cols) in enumerate(self.column_pairs):
            tmp_table = f"tmp_map5_{i}"

            if isinstance(key_cols, str):
                key_expr = f"CAST({key_cols} AS VARCHAR)"
            else:
                key_expr = " || '-' || ".join([f"CAST({c} AS VARCHAR)" for c in key_cols])

            if isinstance(target_cols, str):
                target_expr = f"CAST({target_cols} AS VARCHAR)"
            else:
                target_expr = " || '-' || ".join([f"CAST({c} AS VARCHAR)" for c in target_cols])

            self.con.execute(f"""
                CREATE OR REPLACE TABLE {tmp_table} AS
                SELECT *
                FROM {current_table}
                WHERE {target_expr} IN (
                    SELECT {target_expr}
                    FROM {current_table}
                    GROUP BY {target_expr}
                    HAVING COUNT(DISTINCT {key_expr}) = 1
                )
            """)

            current_table = tmp_table

        self.con.execute(f"CREATE OR REPLACE TABLE map5 AS SELECT * FROM {current_table}")

    def refine_map_from_parquet(self, parquet_path):
        """
        Run the full map refinement pipeline starting from Parquet data.

        Args:
            parquet_path (str): Path to input Parquet file(s).

        Example:
            >>> refiner.refine_map_from_parquet("reads.parquet")
            >>> refiner.refine_map_from_parquet("folder/*.parquet")
        """
        self.create_map1(parquet_path)
        self.create_map2()
        self.create_map3()
        self.create_map4()
        self.create_map5()

    def refine_map_from_db(self):
        """
        Run the map refinement pipeline assuming map1 already exists in the database.

        Example:
            >>> refiner.refine_map_from_db()
        """
        self.create_map2()
        self.create_map3()
        self.create_map4()
        self.create_map5()

    def save_map(self, map_name, map_path):
        """
        Save a DuckDB table to CSV.

        Args:
            map_name (str): Table name to save.
            map_path (str): File path for CSV output.

        Example:
            >>> refiner.save_map("map5", "map5.csv")
        """
        self.con.execute(f"COPY {map_name} TO '{map_path}' WITH (HEADER, DELIMITER ',')")

    def save_loss_table(self, output_csv_path):
        """
        Save a loss table of map1->map5 row counts with descriptive labels and percentages.

        Percentages include:
            - % of previous step
            - % of total reads (map1)

        Args:
            output_csv_path (str): File path to save CSV.

        Example:
            >>> refiner.save_loss_table("map_lengths.csv")
        """
        map_info = [
            ("map1", "Total reads"),
            ("map2", "After removing low quality and undesigned"),
            ("map3", "Grouped counts"),
            ("map4", f"Filtered by reads_threshold > {self.reads_threshold}"),
            ("map5", "Filtered for unique mappings")
        ]

        lengths = []
        existing_tables = [t[0] for t in self.con.execute("SHOW TABLES").fetchall()]

        for table_name, description in map_info:
            if table_name in existing_tables:
                count = self.con.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
            else:
                count = None
            lengths.append({"map": table_name, "description": description, "num_rows": count})

        df = pd.DataFrame(lengths)

        prev_count = None
        map1_count = df.loc[df['map'] == 'map1', 'num_rows'].values[0]
        percent_prev = []
        percent_map1 = []

        for count in df['num_rows']:
            if count is None or prev_count is None:
                percent_prev.append(None)
            else:
                percent_prev.append(round(100 * count / prev_count, 2))
            if count is None or map1_count is None:
                percent_map1.append(None)
            else:
                percent_map1.append(round(100 * count / map1_count, 2))
            prev_count = count if count is not None else prev_count

        df['% of previous step'] = percent_prev
        df['% of total reads'] = percent_map1

        df.to_csv(output_csv_path, index=False)
