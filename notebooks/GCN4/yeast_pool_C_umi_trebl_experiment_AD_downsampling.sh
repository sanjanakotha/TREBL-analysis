#!/bin/bash
#SBATCH --job-name=trebl_downsample
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:05:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/logs/trebl_downsample_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --array=14

# Get the list of FASTQ files
FASTQ_FILES=(
    /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_AD_chunks/*{100,200}_chunks_part_{1..5}.fq.gz
)

# Get the file for this array task
FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

echo "Processing file: $FASTQ_FILE"

# Run TREBL Python script
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/yeast_pool_C_umi_trebl_experiment_AD_downsampling.py "$FASTQ_FILE"
