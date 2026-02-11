#!/bin/bash
#SBATCH --job-name=trebl_downsample
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/logs/trebl_downsample_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --array=0-103
#SBATCH --reservation=maint 


# Load environment
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env

# Get the list of FASTQ files
FASTQ_FILES=(/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_AD_chunks/*.fq.gz)

# Get the file for this array task
FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

echo "Processing file: $FASTQ_FILE"

# Run TREBL Python script
python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/yeast_pool_C_umi_trebl_experiment_AD_downsampling.py "$FASTQ_FILE"
