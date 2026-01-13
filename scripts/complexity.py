import duckdb                # For connecting to your DuckDB database
import pandas as pd          # For DataFrame manipulation
import numpy as np           # For numerical operations (e.g., np.round, np.isfinite)
import seaborn as sns        # For plotting (barplot, styling)
import matplotlib
import os
import matplotlib.pyplot as plt
import tempfile
import os
import shutil
import pathlib
import dask.dataframe as dd

from scripts import preprocess
from scripts import finder

class ComplexityChecker:
    """
    A class to handle sequence preprocessing, barcode extraction, 
    and complexity analysis for sequencing datasets.
    """
    
    def __init__(
        self,
        db_path,
        step_name,
        step1_map_name,
        step_suffix,
        barcode_groups
    ):
        self.con = duckdb.connect(db_path)
        self.step_name = step_name
        self.step1_map_name = step1_map_name
        self.step_suffix = step_suffix
        self.barcode_groups = barcode_groups
        
    def show_tables(self):
        tables = self.con.execute("SHOW TABLES").fetchall()
        for table in tables:
            print(table)
            
    def count_overlap_one_barcode_group(self, barcode_group):
        bc_columns = [barcode.name for barcode in barcode_group]
        table_prefix = self.step_name + "_" + "_".join(bc_columns) + "_"

        # Combine statistics for all barcodes
        map_column_list = ", ".join([f"m.{bc}" for bc in bc_columns])    
        step_column_list = ", ".join([f"s.{bc}" for bc in bc_columns])    
        where_clause = " AND ".join([f"m.{bc} = s.{bc}" for bc in bc_columns])

        pair_query = f"""
            SELECT
                '{",".join(bc_columns)}' AS BC_type,
                (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.step1_map_name} m) t) AS map_unique,
                (SELECT COUNT(*) FROM (SELECT DISTINCT {step_column_list} FROM {table_prefix}{self.step_suffix} s) t1) AS {self.step_name},
                (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.step1_map_name} m JOIN {table_prefix}{self.step_suffix} s ON {where_clause}) t2) AS seen_in_both
        """
        df = self.con.execute(pair_query).df()
        
    
        return df

    def count_overlap(self):
        results = []
        for barcode_group in self.barcode_groups:
            # If not iterable (or is a string), wrap in a list
            if not isinstance(barcode_group, (list, tuple, set)):
                barcode_group = [barcode_group]

            df = self.count_overlap_one_barcode_group(barcode_group)
            results.append(df)
    
        summary_df = pd.concat(results, ignore_index=True)
        summary_df["percent_of_map_seen"] = (
            100 * np.round(summary_df["seen_in_both"] / summary_df["map_unique"], 5)
        )
        return summary_df


