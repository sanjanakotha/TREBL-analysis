import duckdb                # For connecting to your DuckDB database
import pandas as pd          # For DataFrame manipulation
import numpy as np           # For numerical operations (e.g., np.round, np.isfinite)
import seaborn as sns        # For plotting (barplot, styling)
import matplotlib.pyplot as plt  # For figure creation and customization
import preprocess
import finder
import tempfile
import os
import shutil
import pathlib
import dask.dataframe as dd

class ComplexityChecker:
    """
    A class to handle sequence preprocessing, barcode extraction, 
    and complexity analysis for sequencing datasets.
    """
    
    def __init__(
        self,
        db_path,
        bc_objects,
        seq_paths,
        intermediate_output_dir,
        map_name,
        bc_together,
        step_name,
        reverse_complement=True,
        existing_parquet_dir=None
    ):
        """
        Initialize a ComplexityChecker instance.

        Parameters
        ----------
        db_path : str
            Path to existing DuckDB database file with final Step1 map.
        bc_objects : list
            List of barcode objects containing extraction information.
        seq_paths : list
            List of sequencing file paths (TXT/FASTQ/FASTQ.GZ).
        intermediate_output_dir : str
            Directory to store intermediate files.
        map_name : str
            Name of the final Step1 reference mapping table in DuckDB. (ex.) threhsolded or quality_designed)
        bc_together : bool
            Whether to barcodes are in same read or separate.
        step_name : str
            Identifier for the step, used in filenames and table names.
        reverse_complement : bool, default True
            Whether to reverse complement sequences.
        existing_parquet_dir : str or None, default None
            Directory of existing Parquet files to skip preprocessing.
        """
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.seq_paths = seq_paths
        self.bc_objects = bc_objects
        self.reverse_complement = reverse_complement
        self.map_name = map_name.lower().replace(" ", "_")
        self.step_name = step_name.lower().replace(" ", "_")
        self.bc_together = bc_together
        self.bc_columns = [bc_object.name for bc_object in bc_objects]

        # Ensure intermediate output directory exists
        self.intermediate_output_dir = intermediate_output_dir
        os.makedirs(self.intermediate_output_dir, exist_ok=True)

        self.existing_parquet_dir = existing_parquet_dir
            

    def process_step_files(self):
        """
        Load, preprocess, add barcodes, and save sequencing files for this step.

        If `existing_parquet_dir` is provided, loads a Dask DataFrame from Parquet 
        instead of reprocessing raw sequence files.

        Returns
        -------
        df : pd.DataFrame or dask.DataFrame
            Preprocessed dataframe with barcodes added. 
            Returns a Dask DataFrame if loading existing Parquet.
        """
        if not self.existing_parquet_dir:
            # Construct output paths
            parquet_output_dir = os.path.join(
                self.intermediate_output_dir,
                f"{self.step_name}_parquet"
            )
            print(f"Saving {self.step_name} parquet intermediate to {parquet_output_dir}.\n")

            shortened_reads_path = os.path.join(
                self.intermediate_output_dir,
                f"{self.step_name}_reads_shortened.txt"
            )

            # Load and preprocess sequences
            df = preprocess.load_and_shorten_files(
                self.seq_paths, 
                reverse_complement=self.reverse_complement, 
                txt_path=shortened_reads_path
            )

            # Add barcodes
            df = finder.add_multiple_barcodes(self.bc_objects, df)

            # Save as Parquet
            preprocess.save_parquet(df, parquet_output_dir)

        else:
            # Load existing Parquet as Dask DataFrame
            parquet_output_dir = self.existing_parquet_dir
            df = dd.read_parquet(parquet_output_dir)

        # Create DuckDB table for querying
        if self.con is not None:
            self.con.execute(f"""
                CREATE OR REPLACE TABLE {self.step_name}_complexity 
                AS SELECT * FROM read_parquet('{parquet_output_dir}/*')
            """)

        return df

                    
    def complexity_loss_summary(self):
        """
        Summarize barcode complexity losses.

        Returns
        -------
        summary_df : pd.DataFrame
            Summary table containing counts of unique barcodes per column,
            counts passing quality filters, and overlap with reference map.
            Also computes percent of map seen.
        """
        results = []

        # Compute per-barcode statistics
        for bc in self.bc_columns:
            query = f"""
            SELECT
                '{bc}' AS BC_type,
                COUNT(DISTINCT {bc}) AS step1_map,
                (SELECT COUNT(DISTINCT {bc}) FROM {self.step_name}_complexity WHERE {bc} IS NOT NULL) AS {self.step_name},    
                (SELECT COUNT(DISTINCT {bc}) FROM {self.step_name}_complexity WHERE {bc} IS NOT NULL AND {bc}_qual IS TRUE) AS quality_{self.step_name},                
                (SELECT COUNT(DISTINCT s.{bc})
                 FROM {self.step_name}_complexity s
                 JOIN {self.map_name} m
                   ON s.{bc} = m.{bc}) AS seen_in_both
            FROM {self.map_name}
            """
            df = self.con.execute(query).df()
            results.append(df)

        # Combine statistics for all barcodes
        map_column_list = ", ".join([f"m.{bc}" for bc in self.bc_columns])    
        step_column_list = ", ".join([f"s.{bc}" for bc in self.bc_columns])    
        null_conditions = " AND ".join([f"{c} IS NOT NULL" for c in self.bc_columns])
        qual_conditions = " AND ".join([f"{c}_qual IS TRUE" for c in self.bc_columns])
        combined_conditions = f"{null_conditions} AND {qual_conditions}"

        # Compute pairwise overlaps if bc_together is True
        if self.bc_together:
            where_clause = " AND ".join([f"m.{bc} = s.{bc}" for bc in self.bc_columns])
            pair_query = f"""
                SELECT
                    'all_barcodes' AS BC_type,
                    (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.map_name} m) t) AS step1_map,
                    (SELECT COUNT(*) FROM (SELECT DISTINCT {step_column_list} FROM {self.step_name}_complexity s WHERE {null_conditions}) t) AS {self.step_name},
                    (SELECT COUNT(*) FROM (SELECT DISTINCT {step_column_list} FROM {self.step_name}_complexity s WHERE {combined_conditions}) t1) AS quality_{self.step_name},
                    (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.map_name} m JOIN {self.step_name}_complexity s ON {where_clause}) t2) AS seen_in_both
            """
        else:
            where_clause = " AND ".join([f"m.{bc} in (SELECT {bc} FROM {self.step_name}_complexity)" for bc in self.bc_columns])
            pair_query = f"""
            SELECT
                'all_barcodes' AS BC_type,
                (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.map_name} m) t) AS step1_map,
                COUNT(*) AS seen_in_both
            FROM {self.map_name} m
            WHERE {where_clause}
            """

        df_pairs = self.con.execute(pair_query).df()
        results.append(df_pairs)

        # Combine all results and compute percentages
        summary_df = pd.concat(results, ignore_index=True)
        summary_df = summary_df[[
            "BC_type",
            f"step1_map",
            self.step_name,
            "quality_" + self.step_name,
            "seen_in_both",
        ]]
        summary_df["percent_of_map_seen"] = 100 * np.round(summary_df["seen_in_both"] / summary_df[f"step1_map"], 5)

        return summary_df
    
        
    def count_one_bc(self, bc_col_name):
        """
        Count occurrences of a single barcode column in the mapping table,
        limited to barcodes present in the step complexity table.

        Returns
        -------
        df : pd.DataFrame
            Counts of unique values for the barcode column.
        """
        df = self.con.execute(f"""
            SELECT 
                m.{bc_col_name} AS {bc_col_name},
                COUNT(*) AS count
            FROM {self.map_name} AS m
            JOIN {self.step_name}_complexity AS s
                ON m.{bc_col_name} = s.{bc_col_name}
            WHERE m.{bc_col_name} IS NOT NULL
            GROUP BY m.{bc_col_name}
            ORDER BY count DESC
        """).df()
        return df

    def count_all_bc(self):
        """
        Collapse duplicate rows based on all barcode columns and return counts
        per unique combination.

        Returns
        -------
        df_counts : pd.DataFrame
            Counts per unique combination of barcode columns.
        """
        group_cols = ", ".join(self.bc_columns)
        query = f"""
            SELECT {group_cols}, COUNT(*) AS count
            FROM {self.step_name}_complexity
            GROUP BY {group_cols}
            ORDER BY count DESC
        """
        df_counts = self.con.execute(query).df()
        return df_counts.dropna()

    def plot_complexity_loss(self, save_path=None, ax=None, palette="rocket_r"):
        """
        Generate a bar plot visualizing barcode complexity loss.

        Parameters
        ----------
        save_path : str, optional
            File path to save the figure.
        ax : matplotlib.axes.Axes, optional
            Existing axes object to plot on.
        palette : str, optional
            Color palette for seaborn.

        Returns
        -------
        ax : matplotlib.axes.Axes
            Axes object containing the plot.
        """
        plot_df = self.complexity_loss_summary()

        # Melt dataframe for seaborn
        plot_df = plot_df.melt(
            id_vars=["BC_type", "percent_of_map_seen"],
            value_vars=[
                "BC_type",
                f"step1_map",
                self.step_name,
                "quality_" + self.step_name,
                "seen_in_both"
            ],
            var_name="Category",
            value_name="Count"
        )

        # Drop non-numeric counts
        plot_df = plot_df[pd.to_numeric(plot_df["Count"], errors="coerce").notna()]

        # Create figure if not provided
        if ax is None:
            sns.set(style="white", context='talk')
            fig, ax = plt.subplots(figsize=(12, 5), dpi=300)

        # Barplot
        ax = sns.barplot(
            data=plot_df,
            x="Category",
            y="Count",
            hue="BC_type",
            palette=palette
        )

        # Add percent labels for 'seen_in_both'
        seen_in_both_df = plot_df[plot_df["Category"] == "in_seen_in_both"].reset_index()
        for i, row in seen_in_both_df.iterrows():
            height = row["Count"]
            if np.isfinite(height) and np.isfinite(row["percent_of_map_seen"]):
                x = 2 + i / 3.75 - 1 / 3.75
                ax.text(
                    x,
                    height * 0.5,
                    f"{row['percent_of_map_seen']:.1f}%",
                    ha="center",
                    va="bottom",
                    fontsize='xx-small',
                    color="white",
                    fontweight='bold'
                )

        sns.despine()
        plt.xlabel("")
        plt.ylabel("Unique Count")
        plt.legend()
        plt.tight_layout()

        if save_path:
            ax.get_figure().savefig(save_path, bbox_inches='tight')

        return ax


