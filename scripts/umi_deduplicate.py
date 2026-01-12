import duckdb                # For connecting to your DuckDB database
import pandas as pd          # For DataFrame manipulation
import numpy as np           # For numerical operations (e.g., np.round, np.isfinite)
import seaborn as sns        # For plotting (barplot, styling)
import matplotlib.pyplot as plt  # For figure creation and customization
from scripts import preprocess
from scripts import finder
import tempfile
import os
import re
import shutil
import pathlib
import dask.dataframe as dd
from pathlib import Path
import subprocess
os.environ["MPLBACKEND"] = "Agg"
from tqdm import tqdm
from scripts.preprocess import time_it
from scripts import error_correct

class UMIDeduplicator:
    """
    Unified UMI Deduplication class combining simple counts and full UMI-tools pipeline.
    """

    def __init__(self, db_path, bc_objects, step_name, descriptor, step1_map_name,
                 fastq_path, output_path, refined_map_suffix):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.bc_objects = bc_objects
        self.cols = [bc.name for bc in bc_objects]

        self.step_name = step_name
        self.descriptor = descriptor
        self.step1_map_name = step1_map_name
        self.output_path = output_path

        self.refined_map_suffix = refined_map_suffix

        self.table_prefix = self.step_name + "_" + "_".join(self.cols) + "_"
        if descriptor:
            self.table_prefix += f"{descriptor}_"

        # UMI-tools specific attributes
        if fastq_path:
            path = Path(fastq_path)
            base = path.name
            for suffix in [".fastq", ".gz", ".assembled"]:
                base = base.replace(suffix, "")
            self.base = base
            self.fastq_path = fastq_path
 
    def counts_per_umi(self):
        """
        Returns a DataFrame with counts of each UMI per barcode combination.
        Uses the quality_designed table in DuckDB.
        """
        select_cols_sql = ", ".join(self.cols)  # e.g., "ADBC2, HawkBCs"

        query = f"""
            SELECT {select_cols_sql}, UMI, COUNT(*) AS reads
            FROM {self.table_prefix}{self.refined_map_suffix}
            GROUP BY {select_cols_sql}, UMI
            ORDER BY {select_cols_sql}, reads DESC
        """

        df_counts = self.con.execute(query).fetchdf()
        return df_counts
        
    # --------- Methods from SimpleUMIDeduplicator ---------
    def unique_umis_per_barcodes(self):
        """
        Counts unique UMIs per barcode(s), keeping the original columns distinct.
        """
                
        # Include individual columns in SELECT and GROUP BY
        select_cols_sql = ", ".join(self.cols)

        self.new_table_name = f"{self.table_prefix}{self.refined_map_suffix}_umis_collapsed"
        
        #print(f"Counting unique UMIs per barcode(s) and saving as {self.new_table_name}...")
        
        # Build the alias
        #self.alias_name = f"{'_'.join(self.cols)}_umis_unique"
        self.alias_name = 'count'
        
        query = f"""
            CREATE OR REPLACE TABLE {self.new_table_name} AS
            SELECT {select_cols_sql}, 
                   COUNT(DISTINCT UMI) AS {self.alias_name}
            FROM {self.table_prefix}{self.refined_map_suffix}
            GROUP BY {select_cols_sql}
            ORDER BY {self.alias_name} DESC
        """
        
        self.con.execute(query)

    def merge_simple_with_step1_map(self, save = True):
        """
        Performs an inner merge of the current table with step1_map_name
        on the barcode columns (self.cols) and returns the result as a DataFrame.
        """
        select_cols_sql = ", ".join(self.cols)  # e.g., AD_BC, RPTR_BC
        
        join_condition = " AND ".join([f"a.{col} = b.{col}" for col in self.cols])

        query = f"""
            SELECT b.*, a."{self.alias_name}"
            FROM "{self.new_table_name}" AS a
            INNER JOIN "{self.step1_map_name}" AS b
            ON {join_condition}
        """        
        merged_df = self.con.execute(query).fetchdf()
        merged_df["info"] = self.step_name

        filtered_cols = [col for col in merged_df.columns if col != "Designed" and "_qual" not in col and col != "count"]
        merged_df = merged_df[filtered_cols]
            
        if self.output_path and save == True:
            # Ensure the output directory exists
            os.makedirs(self.output_path, exist_ok=True)
        
            # Build a safe filename
            filename = os.path.join(self.output_path, f"{self.base}_{'_'.join(self.cols)}_umis_unique_with_step1_map.csv")
            file_path = os.path.join(self.output_path, filename)
            merged_df.to_csv(file_path, index=False)
            print(f"Saved to {file_path}")
            
        return merged_df

    @time_it
    def save_simple_deduplication(self):
        query = f"""
            SELECT * FROM "{self.new_table_name}"
        """        
        merged_df = self.con.execute(query).fetchdf()

        if self.output_path:
            # Ensure the output directory exists
            os.makedirs(self.output_path, exist_ok=True)

            if self.output_path:
                output_dir = self.output_path
            else:
                output_dir = ""
                
            counts_per_bc = os.path.join(output_dir, f"{self.base}")
            merged_df.to_csv(f"{counts_per_bc}_simple_umi_counts.tsv", index=False, sep = "\t")
            print(f"Saved to {counts_per_bc}_simple_umi_counts.tsv")
            
        return merged_df

    def run_simple_deduplication(self):
        self.unique_umis_per_barcodes()
        return self.merge_simple_with_step1_map()

    def show_tables(self):
        return self.con.execute("SHOW TABLES").fetchdf()

    # --------- Methods from UMIToolsDeduplicator ---------
    @preprocess.time_it
    def generate_fastq(self, suffix="_umi_extracted.fastq"):
        """
        Generate a FASTQ from self.table_prefix.initial:
        - Read sequence = concatenated barcode columns
        - Read ID = UMI
        - Keep all rows (no deduplication)
        - Tracks progress using tqdm
        """
        print("Generating FASTQ with UMIs in header and barcodes as sequence...")
    
        # Concatenate barcode columns
        if len(self.cols) == 1:
            concat_expr = self.cols[0]
        else:
            concat_expr = " || ".join(self.cols)  # no separator
    
        query = f"""
            SELECT UMI, {concat_expr} AS barcode_seq
            FROM {self.table_prefix}{self.refined_map_suffix}
        """
        result = self.con.execute(query).fetchall()
        n_rows = len(result)

        #display(result)
        
        # Prepare output path
        if self.output_path:
            output_file = os.path.join(self.output_path, f"{self.base}{suffix}")
            os.makedirs(self.output_path, exist_ok=True)
        else:
            output_file = f"{self.base}{suffix}"
    
        print(f"Writing FASTQ to {output_file} ({n_rows} reads)...")
    
        with open(output_file, "w") as f:
            for umi, seq in tqdm(result, total=n_rows, desc="Writing FASTQ"):
                # Only keep rows with both a barcode and umi
                if umi is not None and seq is not None and len(seq) > 1 and len(umi) > 1:
                    header = f"@{umi}"
                    plus = "+"
                    qual = "I" * len(seq)  
                    f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
        
        self.umi_fastq = output_file
        print(f"FASTQ complete: {self.umi_fastq}")

    @preprocess.time_it
    def generate_barcode_fasta_and_index(self, suffix = "_barcodes"):
        """
        Save all unique concatenated barcodes as a new table,
        concatenating multiple columns with no separator.
        """
        print("Saving unique barcode(s) as reference file...")
        
        if len(self.cols) == 1:
            concat_expr = self.cols[0]
        else:
            concat_expr = " || ".join(self.cols)  # no separator
    
        new_table_name = f"{self.table_prefix}unique_barcodes"
    
        query = f"""
            CREATE OR REPLACE TABLE {new_table_name} AS
            SELECT DISTINCT {concat_expr} AS barcode
            FROM {self.table_prefix}{self.refined_map_suffix}
        """
    
        print(f"Creating table of unique concatenated barcodes: {new_table_name}")
        self.con.execute(query)

        if self.output_path:
            # Ensure output_path ends with a slash safely
            output_file = os.path.join(self.output_path,f"{self.base}{suffix}.fa")
            

        else:
            output_file = os.path.join(f"{self.base}{suffix}.fa")


        print(f"Writing FASTA to {output_file}...")

        # Fetch all barcodes from DuckDB
        result = self.con.execute(f"SELECT barcode FROM {new_table_name}").fetchall()
        with open(output_file, "w") as f:
            for i, (barcode,) in enumerate(result, start=1):
                # Use the barcode itself as the header
                f.write(f">{barcode}\n{barcode}\n")

        # Index the fasta file too            
        bowtie2_index_prefix = output_file.replace(".fa", "_index")
        self.bowtie2_index_prefix = bowtie2_index_prefix
        
        print(f"Indexing FASTA with bowtie2-build, prefix: {bowtie2_index_prefix}")

        # Run bowtie2-build via module load in a bash shell
        cmd = f"""
        module load bio/bowtie2/2.5.1-gcc-11.4.0
        export BOWTIE2_INDEX_BUILDER_THREADS=32
        bowtie2-build {output_file} {bowtie2_index_prefix}
        """
        
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", stdout=subprocess.DEVNULL)

        os.remove(output_file)

    @preprocess.time_it
    def align_sort_and_deduplicate_umis(self):

        if self.output_path:
            output_dir = self.output_path
            os.makedirs(output_dir, exist_ok=True)
        else:
            output_dir = ""
            
        sam_file = os.path.join(output_dir, f"{self.base}_umi_extracted.sam")
        bam_file = os.path.join(output_dir, f"{self.base}_umi_extracted.bam")
        sorted_bam = os.path.join(output_dir, f"{self.base}_umi_extracted.sorted.bam")
        dedup_bam = os.path.join(output_dir, f"{self.base}_umi_deduplicated.bam")
        #dedup_stats = os.path.join(output_dir, f"{self.base}_deduplicated_stats")
        counts_per_bc = os.path.join(output_dir, f"{self.base}")
        
        # Create the bash script
        script_path = os.path.join(output_dir, f"{self.base}_umi_pipeline.sh")
        
        with open(script_path, "w") as f:
            f.write(f"""#!/bin/bash
    source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/umi_tools
    module load bio/bowtie2/2.5.1-gcc-11.4.0
    module load bio/samtools/1.17-gcc-11.4.0
    
    echo "Aligning .FASTQ to reference .FA ..."
    bowtie2 -p 32 -x {self.bowtie2_index_prefix} -U {self.umi_fastq} -S {sam_file} --norc --no-unal --end-to-end -N 0 -L 12 -k 1 --very-fast
    
    echo "Converting SAM -> BAM ..."
    samtools view -b {sam_file} > {bam_file}
    
    echo "Sorting BAM ..."
    samtools sort -@ 32 {bam_file} -o {sorted_bam}
    
    echo "Indexing BAM ..."
    samtools index {sorted_bam}
    
    echo "Deduplicating UMIs ..."
    umi_tools dedup -I {sorted_bam} -S {dedup_bam}

    echo "Saving final counts..."
    samtools index {dedup_bam}
    umi_tools count -I {dedup_bam} -S {counts_per_bc}_directional_umi_counts.tsv --per-contig &
    wait
    
    echo "Cleaning up intermediate files..."
    rm {sam_file} {bam_file} {sorted_bam} {sorted_bam}.bai {dedup_bam} {dedup_bam}.bai {self.umi_fastq} {self.bowtie2_index_prefix}*

    echo "UMI workflow complete!"
    """)

        # Make the script executable
        os.chmod(script_path, 0o711)
    
        # Run the script
        subprocess.run(f"bash {script_path}", shell=True, check=True, executable="/bin/bash")
    
        os.remove(script_path)

    def run_umi_tools_deduplication(self):
        self.generate_fastq()
        self.generate_barcode_fasta_and_index()
        self.align_sort_and_deduplicate_umis()

    def merge_complex_with_step1_map(self, save = True):
        """
        Merge the UMI-tools deduplicated counts CSV with the step1 map on barcode columns.
        """
        if self.output_path:
            output_dir = self.output_path
        else:
            output_dir = ""
            
        counts_file = os.path.join(output_dir, f"{self.base}_directional_umi_counts.tsv")
        final_counts = pd.read_csv(counts_file, sep="\t")
        
        # Rename columns to match barcode naming convention
        final_counts = final_counts.rename(columns={
            "gene": "_".join(self.cols), 
            "count":  f"{'_'.join(self.cols)}_umis_directional_deduplic"
        })
        
        # Read step1 map from DuckDB
        step1_df = self.con.execute(f'SELECT * FROM "{self.step1_map_name}"').fetchdf()
        
        # Merge on barcode columns
        merge_cols = [col for col in self.cols if col in final_counts.columns and col in step1_df.columns]
        merged_df = step1_df.merge(final_counts, on=merge_cols, how="inner")  
        
        # Optional: filter out unwanted columns like before
        filtered_cols = [col for col in merged_df.columns if col != "Designed" and "_qual" not in col and col != "count"]
        merged_df = merged_df[filtered_cols]
        
        # Add info column
        merged_df["info"] = self.step_name
    
        # Save CSV if output_path specified
        if self.output_path and save == True:
            os.makedirs(self.output_path, exist_ok=True)
            filename = os.path.join(output_dir, f"{self.base}_{'_'.join(self.cols)}_umis_directional_deduplic_with_step1_map.csv")
            file_path = os.path.join(self.output_path, filename)
            merged_df.to_csv(file_path, index=False)
            print(f"Saved to {file_path}")
    
        return merged_df
    
    def run_both_deduplications(self):
        print("Starting simple deduplication.")
        self.unique_umis_per_barcodes()
        print("Finished simple deduplication.\n")
        
        print("Starting UMI Tools directional deduplication.")
        self.run_umi_tools_deduplication()
        print("Finished UMI Tools directional deduplication.\n")

        if self.step1_map_name:
            df1 = self.merge_simple_with_step1_map(save = False)
            df2 = self.merge_complex_with_step1_map(save = False)

            merged_df = pd.merge(df1, df2, how = 'outer')
        
            if self.output_path:
                os.makedirs(self.output_path, exist_ok=True)
                filename = os.path.join(self.output_path, f"{self.base}_{'_'.join(self.cols)}_umis_deduplic_with_step1_map.csv")
                file_path = os.path.join(self.output_path, filename)
                merged_df.to_csv(file_path, index=False)
                print(f"Saved to {file_path}")
                
            return merged_df
        else:
            self.save_simple_deduplication()
            # Need to save the simple deduplication result
            

# Need to test
def run_fastp(input_dir, output_dir, script_path="../savio_jobs/fastp.sh"):
    """
    Run fastp on all .fastq.gz files in input_dir. Outputs as .fastq files.
    Submit the existing fastp.sh script as a SLURM array job.
    
    Args:
        input_dir (str): Path to folder containing .fastq.gz files
        output_dir (str): Path to folder where output should be written
        script_path (str): Path to the existing fastp.sh script
    """
    
    # Resolve absolute paths
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    script_path = os.path.abspath(script_path)
    
    # Check that script exists
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"Script not found: {script_path}")
    
    # Grab all .fastq.gz files to determine array size
    files = sorted([f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")])
    
    if not files:
        raise ValueError(f"No .fastq.gz files found in {input_dir}")
    
    array_range = f"1-{len(files)}"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    
    # Submit the job with SLURM array option
    submit_cmd = [
        "sbatch",
        f"--array={array_range}",
        f"--output={log_dir}/fastp_%A_%a.out",
        script_path,
        input_dir,
        output_dir
    ]
    
    subprocess.run(submit_cmd)
    print(f"Submitted SLURM array job for {len(files)} files using {script_path}.")


def fastp_summary_df(output_dir):
    log_dir = os.path.join(output_dir, "logs")
    # List all .out files
    log_files = [
        f for f in os.listdir(log_dir)
        if f.endswith(".out") and not f.lower().startswith("fastp")
    ]
    
    # Prepare lists
    samples = []
    reads_passed = []
    reads_filtered = []
    
    for f in tqdm(log_files):
        path = os.path.join(log_dir, f)
        with open(path) as fh:
            text = fh.read()
        
        # Sample name (strip _fastp_report.out)
        sample = f.replace("_fastp_report.out","")
        
        # Extract reads passed filter
        match_passed = re.search(r"reads passed filter:\s+([\d,]+)", text)
        match_failed_low = re.search(r"reads failed due to low quality:\s+([\d,]+)", text)
        match_failed_N = re.search(r"reads failed due to too many N:\s+([\d,]+)", text)
        match_failed_short = re.search(r"reads failed due to too short:\s+([\d,]+)", text)
        
        if match_passed:
            passed = int(match_passed.group(1).replace(",",""))
        else:
            passed = 0
        # Total filtered = sum of all failures
        failed = 0
        for m in [match_failed_low, match_failed_N, match_failed_short]:
            if m:
                failed += int(m.group(1).replace(",",""))
        
        samples.append(sample)
        reads_passed.append(passed)
        reads_filtered.append(failed)
    
    # Build dataframe
    df = pd.DataFrame({
        "sample": samples,
        "reads_passed": reads_passed,
        "reads_filtered": reads_filtered
    })
    
    # Sort by total reads if you want
    df["total_reads"] = df["reads_passed"] + df["reads_filtered"]
    df = df.sort_values("total_reads", ascending=False)
    df["filtered_percent"] = df["reads_filtered"] / (df["reads_passed"] + df["reads_filtered"]) * 100
    df = df.reset_index(drop = True)
    return df
