import duckdb
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white", context="talk")
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
import matplotlib.cm as cm

class MapRefiner:
    """
    A class to refine sequencing maps made through BarcodeMapper.

    Args:
        db_path (str): Path for new or existing DuckDB database.
        cols (list[str]): List of barcode or tile column names to process. At least one must be "AD".
        reads_threshold (int): Minimum reads threshold for filtering map5_thresholded.
        column_pairs (list[tuple]): List of (key_cols, target_cols) for enforcing each target maps to one key.
        design_check (bool, optional): Whether to filter out non-designed sequences based on the 'Designed' column.
                                       Defaults to True (keep only Designed sequences).

    Example:
        >>> refiner = MapRefiner(
        ...     db_path="my_database.duckdb",
        ...     cols=["AD", "AD_BC", "RP_BC"],
        ...     reads_threshold=5,
        ...     column_pairs=[("AD", "RP_BC"), (("HawkBCs", "AD_BC"), "RP_BC")],
        ...     design_check=False
        ... )
        >>> refiner.create_initial("reads.parquet")
        >>> refiner.refine_map_from_parquet("reads.parquet")
    """

    # def __init__(self, db_path, cols, reads_threshold, column_pairs, design_check=True):
    #     self.db_path = db_path
    #     self.con = duckdb.connect(self.db_path)
    #     self.cols = cols
    #     self.reads_threshold = reads_threshold
    #     self.column_pairs = column_pairs
    #     self.design_check = design_check

    DEFAULT_MAP_ORDER = [
        "initial",
        "grouped",
        "thresholded",
        "unique_target",
        "quality_designed"
    ]

    def __init__(self, db_path, cols, reads_threshold, column_pairs, design_check=True, map_order=None):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.cols = cols
        self.reads_threshold = reads_threshold
        self.column_pairs = column_pairs
        self.design_check = design_check

        # Use provided map_order or prompt the user
        if map_order is not None:
            self.map_order = map_order
        else:
            # Prompt user for map order
            print("Default map order:")
            for i, name in enumerate(self.DEFAULT_MAP_ORDER, 1):
                print(f"{i}. {name}")

            user_input = input(
                "Enter map order by numbers, comma-separated (e.g., 1,3,2,5,4), or press Enter for default: "
            )

            if user_input.strip() == "":
                self.map_order = self.DEFAULT_MAP_ORDER
            else:
                try:
                    indices = [int(x.strip()) - 1 for x in user_input.split(",")]
                    self.map_order = [self.DEFAULT_MAP_ORDER[i] for i in indices]
                except Exception:
                    print("Invalid input, using default order.")
                    self.map_order = self.DEFAULT_MAP_ORDER

        print("Using map order:", self.map_order)


    def create_initial(self, parquet_path):
        """
        Load raw sequencing data from Parquet file or pattern into initial.

        Args:
            parquet_path (str): Path to input Parquet file(s).

        Example:
            >>> refiner.create_initial("reads.parquet")
            >>> refiner.create_initial("folder/*.parquet")
        """
        self.con.execute(f"""
            CREATE OR REPLACE TABLE initial AS
            SELECT *
            FROM read_parquet('{parquet_path}')
        """)


    def create_quality_designed(self, previous_table = "initial"):
        """
        Remove low-quality or optionally non-designed rows from previous table 
        to create quality_designed.

        Specifically:
            1. Removes low-quality barcodes or tiles (based on *_qual columns).
            2. Optionally removes non-designed sequences if design_check=True.

        Example:
            >>> refiner = MapRefiner(..., design_check=False)
            >>> refiner.create_quality_designed()
        """
        #cols_to_keep = ", ".join(self.cols)
        where_clause = " AND ".join([f"{c}_qual NOT IN (0, FALSE)" for c in self.cols])

        if self.design_check:
            designed_filter = "AND Designed NOT IN (0, FALSE)"
        else:
            designed_filter = ""

        self.con.execute(f"""
            CREATE OR REPLACE TABLE quality_designed AS
            SELECT *
            FROM {previous_table}
            WHERE {where_clause}
            {designed_filter}
        """)

    def create_grouped(self, previous_table="quality_designed"):
        """
        Group by AD + barcode/tile columns, keep count, barcode _qual columns, and Designed.
        """
        # barcode columns are all cols except 'AD'
        qual_cols = [f"{c}_qual" for c in self.cols]

        # select columns: group columns + count + barcode_quals + Designed
        group_cols_sql = ", ".join(self.cols)
        qual_cols_sql = ", ".join(qual_cols + ["Designed"])

        self.con.execute(f"""
            CREATE OR REPLACE TABLE grouped AS
            SELECT {group_cols_sql},
                COUNT(*) AS count,
                {qual_cols_sql}
            FROM {previous_table}
            GROUP BY {group_cols_sql}, {qual_cols_sql}
            ORDER BY count DESC
        """)


    def create_thresholded(self, previous_table = "grouped"):
        """
        Filter grouped by read count threshold to create map5_thresholded.
    
        This now happens BEFORE enforcing unique target constraints.
        """
        self.con.execute(f"""
            CREATE OR REPLACE TABLE thresholded AS
            SELECT *
            FROM {previous_table}
            WHERE count > {self.reads_threshold}
        """)
    
    def create_unique_target(self, previous_table = "thresholded"):
        """
        Enforce mapping constraints from column_pairs on map5_thresholded to create map4_unique_target.
    
        Each target combination should map to exactly one key combination.
        Supports single-column or multi-column keys and targets.
    
        Example:
            >>> refiner.create_map4_unique_target()
        """
    
        for i, (key_cols, target_cols) in enumerate(self.column_pairs):
            tmp_table = f"tmp_map4_{i}"
    
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
                FROM {previous_table}
                WHERE {target_expr} IN (
                    SELECT {target_expr}
                    FROM {previous_table}
                    GROUP BY {target_expr}
                    HAVING COUNT(DISTINCT {key_expr}) = 1
                )
            """)
            previous_table = tmp_table
    
        self.con.execute(f"CREATE OR REPLACE TABLE unique_target AS SELECT * FROM {previous_table}")


    def refine_map_from_parquet(self, parquet_path):
        """
        Run the full map refinement pipeline starting from Parquet data,
        using self.map_order to determine which steps to run.
        """
        # Map table names to functions
        table_to_func = {
            "initial": self.create_initial,
            "quality_designed": self.create_quality_designed,
            "grouped": self.create_grouped,
            "thresholded": self.create_thresholded,
            "unique_target": self.create_unique_target
        }

        prev_table = None

        for table_name in self.map_order:
            func = table_to_func.get(table_name)
            if func is None:
                continue

            # create_initial needs parquet_path, others take previous_table
            if table_name == "initial":
                func(parquet_path)
                prev_table = "initial"
            else:
                func(previous_table=prev_table)
                prev_table = table_name

    def refine_map_from_db(self):
        """
        Run the map refinement pipeline assuming initial already exists in the database,
        using self.map_order to determine which steps to run.
        """
        # Map table names to functions
        table_to_func = {
            "initial": None,  # Already exists
            "quality_designed": self.create_quality_designed,
            "grouped": self.create_grouped,
            "thresholded": self.create_thresholded,
            "unique_target": self.create_unique_target
        }

        # Keep track of the previous table
        prev_table = "initial"

        for table_name in self.map_order:
            func = table_to_func.get(table_name)
            if func is not None:
                func(previous_table=prev_table)
            prev_table = table_name  # Update previous table for next iteration


    def save_map(self, map_name, map_path):
        """
        Save a DuckDB table to CSV.

        Args:
            map_name (str): Table name to save.
            map_path (str): File path for CSV output.

        Example:
            >>> refiner.save_map("map5_thresholded", "map5.csv")
        """
        self.con.execute(f"COPY {map_name} TO '{map_path}' WITH (HEADER, DELIMITER ',')")

    def save_loss_table(self,output_csv_path=False):
        """
        Save a loss table summarizing row counts at each pipeline step (map1 -> map5)
        with descriptive labels and percentages.
    
        Percentages included:
            - % of previous step: fraction of rows remaining compared to the prior map
            - % of total reads (map1): fraction of rows remaining relative to the initial combinations
    
        Args:
            output_csv_path (str or False): File path to save CSV. If False, the table is not saved.
    
        Returns:
            pandas.DataFrame: Table with columns ['map', 'description', 'num_rows', '% of previous step', '% of total reads']
    
        Example:
            >>> df = refiner.save_loss_table("map_lengths.csv")
            >>> print(df)
        """
        map_info_dict = {
            "initial": "Initial combinations",
            "quality_designed": "After removing low quality and undesigned",
            "grouped": "Grouped counts",
            "thresholded": f"Filtered by reads_threshold > {self.reads_threshold}",
            "unique_target": "Filtered for unique targets"
        }

        lengths = []
        existing_tables = [t[0] for t in self.con.execute("SHOW TABLES").fetchall()]

        for table_name in self.map_order:
            if table_name in existing_tables:
                count = self.con.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
            else:
                count = None

            description = map_info_dict[table_name]
            lengths.append({"map": table_name, "description": description, "num_rows": count})

        df = pd.DataFrame(lengths)

        prev_count = None
        map1_count = df.loc[df['map'] == 'initial', 'num_rows'].values[0]
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
        df = df.fillna(100)
        
        if output_csv_path:
            df.to_csv(output_csv_path, index=False)

        return df

    def get_map_df(self, map_name):
        """
        Retrieve any map table as a Pandas DataFrame.
    
        Args:
            map_name (str): Name of the map table to fetch.
                            Should be one of:
                            'initial', 'quality_designed',
                            'grouped', 'map4_unique_target',
                            'map5_thresholded'
    
        Returns:
            pd.DataFrame: DataFrame containing all rows from the specified map.
                          Returns None if the table does not exist.
        """
        existing_tables = [t[0] for t in self.con.execute("SHOW TABLES").fetchall()]
        if map_name not in existing_tables:
            print(f"{map_name} table does not exist.")
            return None
    
        df = self.con.execute(f"SELECT * FROM {map_name}").df()
        return df
        

    def plot_loss(self, save_path=None, ax=None, palette = "rocket_r"):
        """
        Display a horizontal bar plot of the pipeline loss table (map1 -> map5).

        Args:
            save_path (str, optional): File path to save the figure. If None, figure is not saved.
            ax (matplotlib.axes.Axes, optional): Existing matplotlib axis to plot on. 
                                                If None, a new figure and axis are created.

        Returns:
            matplotlib.axes.Axes: The axis containing the plot.
        """
        df = self.save_loss_table()    

        # Set style and create figure if not provided
        if ax is None:
            sns.set(style="white", context='talk')
            fig, ax = plt.subplots(figsize=(6,4), dpi=300)

        # Create  palette for DEFAULT_MAP_ORDER
        cmap = sns.color_palette(palette, n_colors=len(self.DEFAULT_MAP_ORDER))
        palette = {name: cmap[i] for i, name in enumerate(self.DEFAULT_MAP_ORDER)}
        colors = [palette.get(name, (0.5,0.5,0.5)) for name in df['map']]

        # Horizontal bar plot: swap x and y
        sns.barplot(x="num_rows", y="description", data=df, ax=ax, palette=colors)

        # Format y-axis labels
        ax.set_yticklabels(
            df["map"],
        )

        # Add labels at the end of bars
        max_rows = df["num_rows"].max()
        for i, row in df.iterrows():
            x = row["num_rows"]
            ax.text(x, i, f'{row["num_rows"]:,}', va='center', ha='left', fontsize='small')

        # Axis labels
        ax.set_xlabel("Row Count")
        ax.set_ylabel("Map Step")
        
        sns.despine(bottom=True)
        ax.set_xticks([])

        # Optionally save figure
        if save_path:
            ax.get_figure().savefig(save_path, bbox_inches='tight')

        return ax


    def plot_map4_reads(self, save_path=None, ax=None):
        """
        Plot coverage (count distribution) per unique Tile/BC(s) combo in map4.
    
        Args:
            log_x (bool): Whether to use a logarithmic x-axis.
            save_path (str, optional): Path to save the figure. If None, figure is not saved.
            ax (matplotlib.axes.Axes, optional): Existing axis to plot on. If None, a new figure is created.
    
        Returns:
            matplotlib.axes.Axes: The axis containing the plot.
        """
        map4 = self.get_map_df("map4_unique_target")
        
        # Create figure and axes if not provided
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
        
        sns.histplot(map4["count"], bins = 100, log_scale = (True, True), ax = ax)
        sns.despine()

        # Optionally save figure
        if save_path:
            ax.get_figure().savefig(save_path, bbox_inches="tight")
    
        return ax
    
    def _plot_histograms_for_map(self, map_df, axs, color=None, label=None):
        """
        Internal helper: plot read counts and unique BC counts per AD for a given map.
    
        Args:
            map_df (pd.DataFrame): Map dataframe containing at least 'AD', 'count', and barcode columns.
            axs (list of matplotlib.axes.Axes): Axes to plot on. 
                First axis for read counts, subsequent axes for BC counts per AD.
            color (color spec, optional): Color for histogram overlay. Defaults to None.
            label (str, optional): Label for legend. Defaults to None.
    
        Example:
            >>> fig, axs = refiner.plot_map_histograms("grouped")
            >>> refiner._plot_histograms_for_map(map_df, axs, color='blue', label='Map3')
        """
        bc_cols = [c for c in self.cols if c != 'AD']
    
        # 1. Reads per full combination
        sns.histplot(
            map_df["count"],
            bins=100,
            log_scale=(True, True),
            ax=axs[0],
            color=color,
            label=label,
            alpha=0.7,
            edgecolor='none',
            fill=True
        )
    
        # 2. Unique BCs per AD
        for j, bc in enumerate(bc_cols, start=1):
            counts = map_df[['AD', bc]].drop_duplicates().groupby('AD')[bc].size()
            sns.histplot(
                counts,
                bins=100,
                log_scale=(True, True),
                ax=axs[j],
                color=color,
                label=label,
                alpha=0.7,
                edgecolor='none',
                fill=True
            )
    
    
    def plot_map_histograms(self, map_name, axs=None, save=None):
        """
        Plot histograms for a single map, including read counts and unique BC counts per AD.
    
        Args:
            map_name (str): Name of the map dataframe to fetch via `self.get_map_df`.
            axs (matplotlib.axes.Axes or list of Axes, optional): Pre-existing axes to plot on.
                If None, new axes are created.
            save (str, optional): Path to save the figure. If None, figure is not saved.
    
        Returns:
            tuple: (fig, axs)
                fig (matplotlib.figure.Figure or None): The figure object created, or None if axs provided.
                axs (list of matplotlib.axes.Axes): Axes for each subplot.
    
        Example:
            >>> fig, axs = refiner.plot_map_histograms("grouped")
            >>> fig.show()
            >>> refiner.plot_map_histograms("map4_unique_target", save="map4.png")
        """
        map_df = self.get_map_df(map_name)
    
        nplots = 1 + len([c for c in self.cols if c != 'AD'])
        if axs is None:
            fig, axs = plt.subplots(1, nplots, dpi=300, figsize=(5*nplots, 4), sharey=True)
            axs = axs.flatten()
            created_fig = True
        else:
            fig = None
            created_fig = False
    
        # Map description for legend
        map_descriptions = {
            "grouped": "Grouped",
            "map4_unique_target": "Unique targets",
            "map5_thresholded": "Thresholded"
        }
        label = map_descriptions.get(map_name, map_name)
    
        # Plot using helper
        self._plot_histograms_for_map(map_df, axs, color=None, label=label)
    
        # Set axis labels
        axs[0].set_xlabel("+".join(self.cols) + " Reads")
        bc_cols = [c for c in self.cols if c != 'AD']
        for j, bc in enumerate(bc_cols, start=1):
            axs[j].set_xlabel(f"Unique {bc}s per AD")
    
        # Legends and formatting
        axs[-1].legend()
        if created_fig:
            fig.tight_layout(pad=1)
            fig.suptitle(" ".join(map_name.split("_")), y=1.02)
            sns.despine()
            if save is not None:
                fig.savefig(save, dpi=300, bbox_inches="tight")
    
        return fig, axs
    
    
    def plot_maps_3to5_histograms(self, save=False):
        """
        Overlay histograms for maps 3, 4, and 5, showing read counts and unique BC counts per AD.
    
        Args:
            save (bool or str, optional): 
                - False: do not save.
                - True: save to "map3to5_hist.png".
                - str: save to the specified file path.
    
        Returns:
            tuple: (fig, axs)
                fig (matplotlib.figure.Figure): The figure object created.
                axs (list of matplotlib.axes.Axes): Axes for each subplot.
    
        Example:
            >>> refiner.plot_maps_3to5_histograms(save=True) #To display 
            >>> fig, axs = refiner.plot_maps_3to5_histograms(save=True) #To display 
            >>> fig.show()
            >>> refiner.plot_maps_3to5_histograms(save="map_overlay.png") #To save
        """
        maps = {
            "grouped": "Map3",
            "map4_unique_target": "Map4",
            "map5_thresholded": "Map5"
        }
        colors = sns.color_palette("colorblind", n_colors=len(maps))
        map_descriptions = {
            "grouped": "Grouped",
            "map4_unique_target": "Unique targets",
            "map5_thresholded": "Thresholded"
        }
    
        nplots = 1 + len([c for c in self.cols if c != 'AD'])
        fig, axs = plt.subplots(1, nplots, dpi=300, figsize=(5*nplots, 4), sharey=True)
        axs = axs.flatten()
    
        for (map_name, _), color in tqdm(zip(maps.items(), colors), total=len(maps), desc="Processing maps"):
            map_df = self.get_map_df(map_name)
            label = map_descriptions.get(map_name, map_name)
            self._plot_histograms_for_map(map_df, axs, color=color, label=label)
    
        # Set axis labels
        axs[0].set_xlabel("+".join(self.cols) + " Reads")
        bc_cols = [c for c in self.cols if c != 'AD']
        for j, bc in enumerate(bc_cols, start=1):
            axs[j].set_xlabel(f"Unique {bc}s per AD")
    
        # Legends and formatting
        axs[-1].legend()
        fig.tight_layout(pad=1)
        fig.suptitle("Intermediate Maps", y=1.02)
        sns.despine()
    
        if save:
            path = "map3to5_hist.png" if save is True else save
            fig.savefig(path, dpi=300, bbox_inches="tight")
    
        return fig, axs

    def close_connection(self):
        """
        Close the DuckDB connection to release the file lock.

        Example:
            >>> refiner.close_connection()
        """
        if self.con:
            self.con.close()
            self.con = None
