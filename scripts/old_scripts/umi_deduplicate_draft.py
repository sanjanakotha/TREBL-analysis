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
from pathlib import Path
import subprocess
os.environ["MPLBACKEND"] = "Agg"
from tqdm import tqdm

class SimpleUMIDeduplicator:
    """
    UMI Deduplication tools.
    """

    def __init__(self, db_path, bc_objects, step_name, descriptor, step_suffix, step1_map_name, output_path = None):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.cols = [bc_object.name for bc_object in bc_objects]

        # What is the final map from umi step?
        self.step_name = step_name
        self.descriptor = descriptor
        self.step_suffix = step_suffix

        # Get your lookup table
        self.step1_map_name = step1_map_name

        self.table_prefix = self.step_name + "_" + "_".join(self.cols) + "_"
        self.output_path = output_path

    def unique_umis_per_barcodes(self):
        """
        Counts unique UMIs per barcode(s), keeping the original columns distinct.
        """
                
        # Include individual columns in SELECT and GROUP BY
        select_cols_sql = ", ".join(self.cols)

        new_table_name = f"{self.table_prefix}{self.step_suffix}_umis_collapsed"
        
        print(f"Counting unique UMIs per barcode(s) and saving as {new_table_name}...")
        
        # Build the alias
        alias_name = f"unique_{'_'.join(self.cols)}_umis"
        print(alias_name)
        
        query = f"""
            CREATE OR REPLACE TABLE {new_table_name} AS
            SELECT {select_cols_sql}, 
                   COUNT(DISTINCT UMI) AS {alias_name}
            FROM {self.table_prefix}{self.step_suffix}
            GROUP BY {select_cols_sql}
            ORDER BY {alias_name} DESC
        """
        
        self.con.execute(query)
        
    def merge_with_step1_map(self):
        """
        Performs an inner merge of the current table with step1_map_name
        on the barcode columns (self.cols) and returns the result as a DataFrame.
        """
        select_cols_sql = ", ".join(self.cols)  # e.g., AD_BC, RPTR_BC
        
        join_condition = " AND ".join([f"a.{col} = b.{col}" for col in self.cols])

        alias_name = f"unique_{'_'.join(self.cols)}_umis"

        query = f"""
            SELECT b.*, a.{alias_name}
            FROM {self.table_prefix}{self.step_suffix}_umis_collapsed AS a
            INNER JOIN {self.step1_map_name} AS b
            ON {join_condition}
        """
        
        merged_df = self.con.execute(query).fetchdf()
        merged_df["info"] = self.step_name


        if self.output_path:
            # Ensure the output directory exists
            os.makedirs(self.output_path, exist_ok=True)
        
            # Build a safe filename
            filename = f"step1_map_with_{self.step_name}.csv"
            file_path = os.path.join(self.output_path, filename)
        
            # Save CSV

            merged_df.to_csv(file_path, index=False)
            print(f"Saved to {file_path}")
            
        return merged_df

    def get_step1_map_with_unique_umi_counts(self):
        self.unique_umis_per_barcodes()
        return self.merge_with_step1_map()

class UMIToolsDeduplicator:

    def __init__(self, db_path, bc_objects, step_name, descriptor, step1_map_name, fastq_path, umi_length, output_path = None):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.cols = [bc_object.name for bc_object in bc_objects]
        self.bc_objects = bc_objects
        
        # What is the final map from umi step?
        self.step_name = step_name
        self.descriptor = descriptor

        # Get your lookup table
        self.step1_map_name = step1_map_name

        self.table_prefix = self.step_name + "_" + "_".join(self.cols) + "_"
        
        if descriptor:
            self.table_prefix = f"{self.table_prefix}{descriptor}_"
            
        path = Path(fastq_path)
        
        # Strip known suffixes
        base = path.name
        for suffix in [".fastq", ".gz", ".assembled"]:
            base = base.replace(suffix, "")
        self.base = base

        self.output_path = output_path
        
        self.fastq_path = fastq_path
        self.umi_length = umi_length

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
            FROM {self.table_prefix}quality_designed
        """
        result = self.con.execute(query).fetchall()
        n_rows = len(result)

        #display(result)
        
        # Prepare output path
        if self.output_path:
            output_file = os.path.join(self.output_path, f"{self.base}{suffix}")
        else:
            output_file = f"{self.base}{suffix}"
    
        print(f"Writing FASTQ to {output_file} ({n_rows} reads)...")
    
        with open(output_file, "w") as f:
            for umi, seq in tqdm(result, total=n_rows, desc="Writing FASTQ"):
                # Only keep rows with both a barcode and umi
                if umi is not None and seq is not None and len(seq) > 1 and len(umi) > 1:
                    header = f"@{umi}"
                    plus = "+"
                    qual = "I" * len(seq)  # dummy quality
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
            FROM {self.table_prefix}quality_designed
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

    def run_full_umi_pipeline(self):
        self.generate_fastq()
        self.generate_barcode_fasta_and_index()
        self.align_sort_and_deduplicate_umis()

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
    
    # Submit the job with SLURM array option
    submit_cmd = [
        "sbatch",
        f"--array={array_range}",
        script_path,
        input_dir,
        output_dir
    ]
    
    subprocess.run(submit_cmd)
    print(f"Submitted SLURM array job for {len(files)} files using {script_path}.")
