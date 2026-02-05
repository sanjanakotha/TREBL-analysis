import matplotlib.pyplot as plt
from pathlib import Path
import duckdb
import os
import pandas as pd
from scripts import initial_map, map_refiner, complexity, finder, preprocess, error_correct, plotting, umi_deduplicate

import math
import seaborn as sns
import re
from tqdm import tqdm
import re


# Defines the sequence of refinement operations for each step.
# Indexed by step name and whether error correction is enabled.
MAP_ORDERS = {
    "step1": {
        False: ['grouped','thresholded','barcode_exists','unique_target','quality','designed'],
        True: ['barcode_exists', 'quality', 'error_corrected','grouped','thresholded','unique_target','designed'],
    },
    "step2": {
        False: ['grouped','thresholded','barcode_exists', 'quality','designed'],
        True: ['barcode_exists', 'quality','error_corrected','grouped','thresholded','designed'],
    },
}


class TreblPipeline:
    """
    End-to-end TREBL analysis steps.

    Attributes:
        db_path (str): Path to the DuckDB database used for all steps.
        error_correction (bool): Whether to apply error-correction
            operations during refinement.
        output_figures_path (str | None): Directory to write figures.
            If None, figures are not generated.
    """

    def __init__(
        self,
        db_path,
        design_file_path,
        error_correction=False,
        output_path=None,
        test_n_reads=False
    ):
        """
        Initialize a TREBL pipeline run.

        Args:
            db_path (str): Path to the DuckDB database.
            design_file_path (str | None): Path to the design file. If None,
                design-based filtering is skipped.
            error_correction (bool, optional): Whether to enable error correction
                during refinement. Defaults to False.
            output_path (str | None, optional): Output directory for results
                and figures. If None, nothing is written to disk.
            test_n_reads (int | bool, optional): If set, limits mapping to the
                first N reads for testing/debugging.
        """
        self.con = duckdb.connect(db_path)
        self.db_path = db_path
        self.error_correction = error_correction
        self.design_file_path = design_file_path
        self.test_n_reads = test_n_reads
        
        if output_path:
            self.output_path = Path(output_path)
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_figures_path = self.output_path / "figures"
            self.output_figures_path.mkdir(exist_ok=True)
        else:
            self.output_path = None
            self.output_figures_path = None
            
    def _run_initial_mappers(self, mapper_specs):
        """
        Helper function to run one or more InitialMapper instances.

        Args:
            mapper_specs (list[dict]): List of mapper specifications. Each dict
                must contain:
                    - seq_file (str): FASTQ file path
                    - step_name (str): Step name prefix for DuckDB tables
                    - bc_objects (list): Barcode objects
                    - reverse_complement (bool): Whether to reverse complement reads
                    - design_file_path (str | None): Design file path
        """
        ...
        for spec in mapper_specs:
            mapper = initial_map.InitialMapper(
                db_path=self.db_path,
                seq_file=spec["seq_file"],
                step_name=spec["step_name"],
                bc_objects=spec["bc_objects"],
                reverse_complement=spec["reverse_complement"],
                design_file_path=spec["design_file_path"],
            )
            if self.test_n_reads:
                mapper.create_test_map(test_n_reads=self.test_n_reads)
            else:
                mapper.create_map()

    def _run_refiners(self, refiners, plot_titles):
        """
        Helper function to run MapRefiner instances and optionally generate loss plots.

        Args:
            refiners (list[MapRefiner]): Refiners to execute.
            plot_titles (list[str | None]): Titles for loss plots.

        Returns:
            list[MapRefiner]: Executed refiners.
        """
        for refiner, plot_title in zip(refiners, plot_titles):
            refiner.refine_map_from_db()
    
            if self.output_figures_path and plot_titles:
                refiner.plot_loss(text_offset = -0.2)
                if plot_title:
                    plt.title(plot_title)
    
        return refiners

    def reads_distribution(self, seq_file, bc_objects, step_name, reverse_complement):
        """
        Visualize read-count distributions for a single FASTQ file to help decide a threshold. 

        This performs initial mapping and grouping and plots a reads histogram.

        Args:
            seq_file (str): FASTQ file path.
            bc_objects (list): Barcode objects.
            step_name (str): Step name used for DuckDB tables.
            reverse_complement (bool): Whether to reverse complement reads.
        """
        
        self._run_initial_mappers([
            {
                "seq_file": seq_file,
                "step_name": step_name,
                "bc_objects": bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": self.design_file_path,
            }
        ])

        # Create refiner
        refiner = map_refiner.MapRefiner(
            db_path=self.db_path,
            bc_objects=bc_objects,
            column_pairs=[],
            reads_threshold=1,
            map_order=["grouped"],
            step_name=step_name,
            output_figures_path=self.output_figures_path,
            design_file = self.design_file_path
        )

        refiners = self._run_refiners([refiner], plot_titles=[f"{step_name} grouped"])
        
    def step1_reads_distribution(self, seq_file, bc_objects, reverse_complement, step_suffix=""):
        """
        Convenience wrapper for Step 1 read distributions.

        Args:
            seq_file (str): FASTQ file path.
            bc_objects (list): Barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            step_suffix (str, optional): Suffix appended to the step name.
        """
        step_name = "step1" + step_suffix
        self.reads_distribution(seq_file, bc_objects, step_name, reverse_complement)

    def step2_reads_distribution(self, 
                                 AD_seq_file,
        AD_bc_objects,
        RT_seq_file,
        RT_bc_objects,
        reverse_complement,
        step_suffix=""):
        """
        Convenience wrapper for Step 2 AD and RT libraries.

        Args:
            AD_seq_file (str): FASTQ file for AD reads.
            AD_bc_objects (list): AD barcode objects.
            RT_seq_file (str): FASTQ file for RT reads.
            RT_bc_objects (list): RT barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            step_suffix (str, optional): Suffix appended to the step name.
        """

        step_name = "step2" + step_suffix
        
        self.reads_distribution(AD_seq_file, AD_bc_objects, step_name, reverse_complement)
        self.reads_distribution(RT_seq_file, RT_bc_objects, step_name, reverse_complement)

    def trebl_experiment_reads_distribution(self, AD_seq_files,
        AD_bc_objects,
        RT_seq_files,
        RT_bc_objects,
        reverse_complement,
        step_suffix=""):
        """
        Convenience wrapper for TREBL Experiment AD and RT libraries.

        Args:
            AD_seq_file (str): FASTQ file for AD reads.
            AD_bc_objects (list): AD barcode objects.
            RT_seq_file (str): FASTQ file for RT reads.
            RT_bc_objects (list): RT barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            step_suffix (str, optional): Suffix appended to the step name.
        """

        step_name_base = "trebl_experiment" + step_suffix
        
        for file_path in AD_seq_files:
            base_name = os.path.basename(file_path)
            name_only = base_name.split('.')[0]
            step_name = f"{step_name_base}_{name_only}"
            self.reads_distribution(file_path, AD_bc_objects, step_name, reverse_complement)
            
        for file_path in RT_seq_files:
            base_name = os.path.basename(file_path)
            name_only = base_name.split('.')[0]
            step_name = f"{step_name_base}_{name_only}"
            self.reads_distribution(file_path, RT_bc_objects, step_name, reverse_complement)
            
    def run_step_1(
        self,
        seq_file,
        bc_objects,
        column_pairs,
        reads_threshold,
        reverse_complement=False,
        step_suffix = ""
    ):
        """
        Run TREBL Step 1: map reads to designed AD barcodes.

        Performs initial mapping, refinement, optional error correction,
        and outputs the final designed map.

        Args:
            seq_file (str): FASTQ file containing Step 1 reads.
            bc_objects (list): Barcode objects defining extraction rules.
            column_pairs (list[tuple]): Column-pair checks used
                to remove barcode collisions. 

                Each tuple has the form:

                    (key_column(s), target_column(s))

                where `key_column(s)` and `target_column(s)` can be either
                a single column name (str) or a tuple/list of column names.

                The constraint enforces that, for each **key**, at least
                90% of its reads must map to a single **target**. If a key
                value maps to multiple targets without one reaching the 90%
                threshold, the key is considered ambiguous and discarded.

                **Single-column example**:

                    column_pairs = [("RPTR_BC", "AD")]

                This ensures that each reporter barcode (``RPTR_BC``) maps
                predominantly to a single AD. If a reporter barcode has
                reads mapping to multiple ADs without one exceeding 90%,
                that reporter barcode is removed.

                **Multi-column example**:

                    column_pairs = [
                        (("RPTR_BC"), ("Hawk_BC", "AD_BC"))
                    ]

                This checks that each RPTR BC
                maps to a single Hawkins AD barcode and AD barcode combination 
                with ≥90% of reads. Ambiguous key combinations are removed.

                Multiple constraints may be applied sequentially:

                    column_pairs = [
                        ("AD_BC", "AD"),
                        ("RPTR_BC", "AD"),
                        (("AD_BC", "RPTR_BC"), ("AD", "SampleID")),
                    ]

            reads_threshold (int): Minimum number of reads required for a
                barcode or barcode pair to be kept.
            reverse_complement (bool, optional): Whether reads should be
                reverse complemented prior to barcode extraction.
                Defaults to False.
            step_suffix (str, optional): Suffix appended to the step name
                for DuckDB table naming. This is useful when you want
                to distinguish multiple runs or subsets of the same step.

                Example:
                    "_spike_in" — if processing spike-in samples separately
                    from your main dataset.
                    "_new_data" — for data from new sequencing run.

        Returns:
            pd.DataFrame: Final designed Step 1 mapping after all refinement
            steps have been applied.
        """


        step_name = "step1" + step_suffix

        
        # Initial mapping
        self._run_initial_mappers([
            {
                "seq_file": seq_file,
                "step_name": step_name,
                "bc_objects": bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": self.design_file_path,
            }
        ])

        manual_ec_threshold = True
        if self.error_correction and reads_threshold == 1:
                manual_ec_threshold = False

        # Create refiner
        refiner = map_refiner.MapRefiner(
            db_path=self.db_path,
            bc_objects=bc_objects,
            column_pairs=column_pairs,
            reads_threshold=reads_threshold,
            map_order=MAP_ORDERS["step1"][self.error_correction],
            step_name=step_name,
            output_figures_path=self.output_figures_path,
            manual_ec_threshold = manual_ec_threshold,
            design_file = self.design_file_path
        )
    
        # Use _run_refiners to handle refinement and plotting
        refiners = self._run_refiners([refiner], plot_titles=["Step 1"])
    
        # Get the dataframe
        df = refiners[0].get_map_df('designed')
        if self.output_path:
            df.to_csv(self.output_path / f"{step_name}.csv", index=False)

        self.step1_df = df
        
        return df

    def run_step_2(
        self,
        AD_seq_file,
        AD_bc_objects,
        RT_seq_file,
        RT_bc_objects,
        reverse_complement,
        reads_threshold_AD,
        reads_threshold_RT,
        step1_map_name,
        step_suffix="",
    ):
        """
        Run TREBL Step 2: Analyze intermediate complexity.

        Performs separate refinement for AD and RT libraries and computes
        overlap with the Step 1 map.

        Args:
            AD_seq_file (str): FASTQ file containing AD reads.
            AD_bc_objects (list): Barcode objects for ADs.
            RT_seq_file (str): FASTQ file containing RT reads.
            RT_bc_objects (list): Barcode objects for reporters.
            reverse_complement (bool): Whether to reverse complement reads.
            reads_threshold_AD (int): Minimum reads per AD barcode.
            reads_threshold_RT (int): Minimum reads per RT barcode.
            step1_map_name (str): Name of the Step 1 DuckDB table.
            step_suffix (str, optional): Suffix appended to the step name.

        Returns:
            dict: Dictionary with keys:
                - "AD_step2": AD DataFrame
                - "RT_step2": RT DataFrame
                - "step1_overlap": Overlap statistics
        """
        step_name = "step2" + step_suffix
        
        # --- Initial mapping ---
        self._run_initial_mappers([
            {
                "seq_file": AD_seq_file,
                "step_name": step_name,
                "bc_objects": AD_bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": self.design_file_path,  # AD uses design
            },
            {
                "seq_file": RT_seq_file,
                "step_name": step_name,
                "bc_objects": RT_bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": None,  # RT skips design
            },
        ])
    
        refiners = []
        
        for bc_objs, reads_threshold in zip(
            (AD_bc_objects, RT_bc_objects),
            (reads_threshold_AD, reads_threshold_RT)
        ):
            step2_map_order = MAP_ORDERS["step2"][self.error_correction].copy()
            
            refiners.append(
                map_refiner.MapRefiner(
                    db_path=self.db_path,
                    bc_objects=bc_objs,
                    column_pairs=[],
                    reads_threshold=reads_threshold,
                    map_order=step2_map_order,
                    step_name=step_name,
                    design_file=self.design_file_path if bc_objs is AD_bc_objects else None,
                    output_figures_path=self.output_figures_path,
                )
            )
            
        # --- Run refiners ---
        refiners = self._run_refiners(
            refiners,
            plot_titles=["AD Step 2", "RT Step 2"],
        )
    
        # --- Compute AD–reporter complexity and overlap ---
        checker = complexity.ComplexityChecker(
            db_path=self.db_path,
            step_name=step_name,
            step1_map_name=step1_map_name,
            step_suffix="designed",
            barcode_groups=[AD_bc_objects, RT_bc_objects],
        )
    
        overlap = checker.count_overlap()
    
        AD_df = refiners[0].get_map_df('designed')
        RT_df = refiners[1].get_map_df('designed')
    
        # --- Save CSVs ---
        if self.output_path:
            AD_df.to_csv(self.output_path / f"{step_name}_AD.csv", index=False)
            RT_df.to_csv(self.output_path / f"{step_name}_RT.csv", index=False)
    
        return {
            "AD_step2": AD_df,
            "RT_step2": RT_df,
            "step1_overlap": overlap
        }

    
    def _duckdb_safe_name(self, base_name):
        name = base_name.replace("-", "_").replace(" ", "_")
        name = re.sub(r'[^0-9a-zA-Z_]', '_', name)
        if re.match(r'^\d', name):
            name = f"_{name}"
        return name
    

    def _run_trebl_experiment_helper(
        self,
        seq_files,
        bc_objects,
        reverse_complement,
        reads_threshold=1,
        umi_object=None,
        step_name_suffix="",
        count_col_name=None,
        gene_col_name=None,
        concat_gene=False,
    ):
        """
        Core TREBL experiment runner.

        Automatically selects UMI or non-UMI workflows based on whether
        a UMI object is provided.

        Args:
            seq_files (list[str]): FASTQ file paths.
            bc_objects (list): Barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            reads_threshold (int, optional): Minimum reads threshold.
            umi_object (optional): UMI configuration object.
            step_name_suffix (str, optional): Suffix appended to step names.
            count_col_name (str, optional): Column name for UMI counts.
            gene_col_name (str, optional): Column name for gene identifiers.
            concat_gene (bool, optional): Whether to concatenate barcode
                columns to form a gene identifier.

        Returns:
            pd.DataFrame: Aggregated experiment results.
        """
        
        step_name_prefix = "trebl_experiment_" + step_name_suffix
        
        results = []
        simple_results = []
    
        for file_path in seq_files:
            base_name = os.path.basename(file_path)
            name_only = base_name.split('.')[0]
            name_only = self._duckdb_safe_name(name_only)
            step_name = f"{step_name_prefix}{name_only}"
    
            output_dir = self.output_path / step_name
            output_dir.mkdir(parents=True, exist_ok=True)
    
            design_file_path = self.design_file_path if "AD" in [bc.name for bc in bc_objects] else None
    
            # --- Initial mapping ---
            mapper_kwargs = dict(
                db_path=self.db_path,
                step_name=step_name,
                seq_file=file_path,
                bc_objects=bc_objects,
                reverse_complement=reverse_complement,
                design_file_path=design_file_path
            )
            if umi_object:
                mapper_kwargs["umi_object"] = umi_object
    
            mapper = initial_map.InitialMapper(**mapper_kwargs)
            mapper.create_map()
    
            # --- Refinement ---
            manual_ec_threshold = not (self.error_correction and reads_threshold == 1)
            map_order = ["quality", "error_corrected"] if self.error_correction else ["quality"]
            if not umi_object:
                # non-UMI workflow adds additional steps
                map_order = map_order + (["grouped", "thresholded", "designed"] if not self.error_correction else ["grouped", "thresholded", "designed"])
    
            refiner = map_refiner.MapRefiner(
                db_path=self.db_path,
                bc_objects=bc_objects,
                column_pairs=[],
                reads_threshold=reads_threshold,
                map_order=map_order,
                step_name=step_name,
                output_figures_path=output_dir,
                manual_ec_threshold=manual_ec_threshold,
            )
            refiner.refine_map_from_db()
            refiner.plot_loss()
            if self.error_correction:
                refiner.plot_error_correction()
    
            if umi_object:
                # --- Deduplication ---
                refined_map_suffix = "error_corrected" if self.error_correction else "quality"
                deduplicator = umi_deduplicate.UMIDeduplicator(
                    db_path=self.db_path,
                    bc_objects=bc_objects,
                    step_name=step_name,
                    descriptor="",
                    step1_map_name=None,
                    fastq_path=file_path,
                    output_path=output_dir,
                    refined_map_suffix=refined_map_suffix,
                )
                deduplicator.run_both_deduplications()
    
                # --- Load results ---
                complex_df = pd.read_csv(output_dir / f"{name_only}_directional_umi_counts.tsv", sep="\t")
                simple_df = pd.read_csv(output_dir / f"{name_only}_simple_umi_counts.tsv", sep="\t")
                complex_df["name"] = name_only
                simple_df["name"] = name_only
                results.append(complex_df)
                simple_results.append(simple_df)

                # Reads per UMI
                one_file_reads_per_UMI = deduplicator.counts_per_umi()
                one_file_reads_per_UMI["name"] = name_only
                one_file_reads_per_UMI.to_csv(output_dir / f"{name_only}_reads_per_umi.tsv", sep="\t")

        # --- Merge results ---
        if umi_object:
            complex_df = pd.concat(results, ignore_index=True).rename(columns={"count": count_col_name, "gene": gene_col_name})
            simple_df = pd.concat(simple_results, ignore_index=True).rename(columns={"count": count_col_name})
    
            if concat_gene:
                concat_cols = [bc.name for bc in bc_objects]
                simple_df[gene_col_name] = simple_df[concat_cols].agg("".join, axis=1)
    
            merged = pd.merge(
                complex_df,
                simple_df,
                on=[gene_col_name, "name"],
                suffixes=("_complex", "_simple"),
            )
            return merged
        else:
             # --- Non-UMI workflow: grab all relevant tables ---
            tables = refiner.show_tables()
            first_bc_name = bc_objects[0].name
            for table in tables:
                if step_name_prefix in table[0] and first_bc_name in table[0]:
                    df = refiner.get_map_df(table[0])
                    df["sample"] = table[0][len(step_name_prefix):]
                    results.append(df)

            return pd.concat(results, ignore_index=True)

    def trebl_experiment_analysis(
        self,
        AD_seq_files,
        AD_bc_objects,
        RT_seq_files,
        RT_bc_objects,
        reverse_complement,
        step1_map_name=None,
        AD_umi_object=None,
        RT_umi_object=None,
        reads_threshold_AD=1,
        reads_threshold_RT=1,
        step_name_suffix="",
    ):
        """
        Run TREBL experiment analysis for both AD and RT libraries.

        The workflow automatically selects between a UMI or non-UMI
        pipeline depending on whether a UMI object is provided.

        Args:
            AD_seq_files (list[str]): Paths to FASTQ files containing AD reads.
            AD_bc_objects (list): Barcode objects for AD library extraction.
            RT_seq_files (list[str]): Paths to FASTQ files containing reporter reads.
            RT_bc_objects (list): Barcode objects for reporter library extraction.
            reverse_complement (bool): Whether reads should be reverse complemented
                prior to barcode extraction.
            step1_map_name (str): DuckDB table name from Step 1 mapping used to
                compute overlap.
            AD_umi_object (optional): UMI object for AD library. If provided,
                UMI deduplication is performed.
            RT_umi_object (optional): UMI object for RT library. If provided,
                UMI deduplication is performed.
            reads_threshold_AD (int, optional): Minimum reads per AD barcode
                to be retained. Defaults to 1.
            reads_threshold_RT (int, optional): Minimum reads per RT barcode
                to be retained. Defaults to 1.
            step_name_suffix (str, optional): Suffix for DuckDB tables and
                output naming after "trebl_experiment_".

        Returns:
            dict: Dictionary containing the final merged results for AD and RT.
                Keys:
                    - "AD_results" (pd.DataFrame): TREBL experiment results for AD
                        library, with barcodes and reads or unique UMI counts. 
                    - "RT_results" (pd.DataFrame): TREBL experiment results for reporter
                        library, with barcodes and reads or unique UMI counts. 

        Notes:
            - If UMI objects are provided, UMI-based deduplication is applied
              and both "complex" and "simple" UMI count tables are merged.
            - If `output_path` is set, results are saved as CSV files:
                "{output_path}/AD_trebl_experiment_results.csv" and
                "{output_path}/RT_trebl_experiment_results.csv".
            - The method also generates barcode quality/loss plots/
        """
        
        step_name_prefix = "trebl_experiment_" + step_name_suffix
        
        experiments = {
            "AD": {
                "seq_files": AD_seq_files,
                "bc_objects": AD_bc_objects,
                "umi_object": AD_umi_object,
                "reads_threshold": reads_threshold_AD,
                "count_col_name": "AD_umi_count",
                "gene_col_name": "AD_ADBC_concat",
                "concat_gene": True,
                "output_file": "ADBC_trebl_experiment_results.csv",
            },
            "RT": {
                "seq_files": RT_seq_files,
                "bc_objects": RT_bc_objects,
                "umi_object": RT_umi_object,
                "reads_threshold": reads_threshold_RT,
                "count_col_name": "RTBC_umi_count",
                "gene_col_name": "RTBC",
                "concat_gene": False,
                "output_file": "RTBC_trebl_experiment_results.csv",
            },
        }
    
        results = {}
    
        for name, spec in experiments.items():
            df = self._run_trebl_experiment_helper(
                seq_files=spec["seq_files"],
                bc_objects=spec["bc_objects"],
                reverse_complement=reverse_complement,
                reads_threshold=spec["reads_threshold"],
                umi_object=spec["umi_object"],
                count_col_name=spec.get("count_col_name"),
                gene_col_name=spec.get("gene_col_name"),
                concat_gene=spec.get("concat_gene", False),
                step_name_suffix=step_name_suffix
            )
    
            if self.output_path:
                df.to_csv(self.output_path / f"{name}_trebl_experiment_results.csv", index=False)
    
            results[f"{name}_results"] = df

        if step1_map_name:            
            self.plot_trebl_experiment_loss(
                AD_bc_objects, 
                step1_map_name, 
                step_name_prefix=step_name_prefix)
    
            self.plot_trebl_experiment_loss(
                RT_bc_objects, 
                step1_map_name, 
                step_name_prefix=step_name_prefix)
        
        return results

    def plot_trebl_experiment_loss(self, bc_objects, step1_map_name=None, step_name_prefix="trebl_experiment_"):
        """
        Plot barcode quality and mapping loss for a TREBL experiment.

        Generates bar plots showing:
            1. Total number of reads in each initial mapping table.
            2. Number of reads passing barcode quality checks (`_qual` columns).
            3. Number of reads that match the Step 1 map for all barcodes.

        Plots are saved as a PNG file in `self.output_figures_path` and 
        returned as Matplotlib figure and axes objects.

        Args:
            bc_objects (list): List of barcode objects to evaluate.
            step1_map_name (str): Name of the Step 1 DuckDB table used to
                calculate overlap with mapped barcodes.
            step_name_prefix (str, optional): Prefix used to identify TREBL
                experiment tables in DuckDB. Defaults to "trebl_experiment_".

        Returns:
            tuple: (fig, axes)
                - fig (matplotlib.figure.Figure): Figure object containing all subplots.
                - axes (list[matplotlib.axes._subplots.AxesSubplot]): Flattened list of subplot axes.
                  Returns None, None if no matching tables are found.

        """
        
        # Connect to DuckDB
        con = duckdb.connect(self.db_path)
    
        # Get all tables matching step_name_prefix and bc_object names
        tables = con.execute("SHOW TABLES").fetchall()
        bc_names = [bc.name for bc in bc_objects]
    
        result_prefixes = [
            table[0] for table in tables
            if step_name_prefix in table[0] 
            and any(bc_name in table[0] for bc_name in bc_names)
            and 'initial' in table[0]
        ]
        
                
        if not result_prefixes:
            print("No matching tables found.")
            return None, None
    
        num_plots = len(result_prefixes)
        max_cols = 5
        ncols = min(max_cols, num_plots)
        nrows = math.ceil(num_plots / ncols)
    
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3), dpi=300, sharey=True, sharex=True)
        axes = axes.flatten() if num_plots > 1 else [axes]
        
        bc_names_str = "_".join([bc.name for bc in bc_objects])

        for i, file_name in tqdm(enumerate(result_prefixes), total=len(result_prefixes), desc="Plotting BCs"):
            # Total count
            total_count = con.execute(f'SELECT COUNT(*) FROM "{file_name}"').fetchone()[0]
    
            # Count rows where all BCs are qualified (if the _qual column exists)
            qual_cols = [f'"{bc.name}_qual"' for bc in bc_objects]

            if qual_cols:
                qual_conditions = " AND ".join(qual_cols)
                both_true = con.execute(f'SELECT COUNT(*) FROM "{file_name}" WHERE {qual_conditions}').fetchone()[0]
            else:
                both_true = 0
    
            # Count rows present in Step1 map for all BCs
            join_conditions = " AND ".join([f'm."{bc.name}" = s."{bc.name}"' for bc in bc_objects])

            if step1_map_name== None:
                step1_count = 0
            else:
                step1_count = con.execute(f'''
                    SELECT COUNT(*) 
                    FROM "{file_name}" AS m
                    JOIN "{step1_map_name}" AS s
                    ON {join_conditions}
                ''').fetchone()[0]
        
            # Prepare DataFrame for plotting
            plot_counts = pd.DataFrame({
                "Category": ["Total", "BC\nQual", "In\nStep1"],
                "Count": [total_count, both_true, step1_count]
            })
    
            ax = axes[i]
            sns.barplot(x="Category", y="Count", data=plot_counts, ax=ax, palette=["gray", "green", "blue"])
            for container in ax.containers:
                ax.bar_label(container, fmt='%d', label_type='edge', fontsize='small', padding=2)

            title_regex = f'{step_name_prefix}(.*)_{bc_names_str}_initial'
            match = re.search(title_regex, file_name)
            group_name = match.group(1) if match else file_name
            ax.set_title(str(group_name), fontsize='medium', y = 1.1)
            ax.set_xlabel("")
            ax.set_ylabel("")
    
        # Hide unused axes
        for j in range(num_plots, len(axes)):
            axes[j].set_visible(False)
    
        fig.supylabel("Count")
        sns.despine()
        plt.tight_layout(pad=1)

        fig.suptitle("Trebl Experiment Loss")
    
        save_path = os.path.join(self.output_figures_path, f"{step_name_prefix}{bc_names_str}_bc_quality_loss.png")
        plt.savefig(save_path)

        con.close()

    
        return fig, axes


