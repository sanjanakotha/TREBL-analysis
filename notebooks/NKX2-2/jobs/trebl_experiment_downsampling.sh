#!/bin/bash
#SBATCH --job-name=trebl_downsample
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/NKX2-2/logs/trebl_downsample_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --array=0

# -------------------------
# Load conda and environment
# -------------------------
# module load anaconda3/2024.02-1-11.4
# source $(conda info --base)/etc/profile.d/conda.sh  # initialize conda
# conda activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env
# set -e  # stop the job if any command fails

# -------------------------
# Create array of FASTQ files
# -------------------------
FASTQ_FILES=(/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/NKX2-2_trebl_exp_chunks/AD_puro_only/*.fq.gz \
             /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/NKX2-2_trebl_exp_chunks/RT_puro_only/*.fq.gz)

echo "Found ${#FASTQ_FILES[@]} FASTQ files"
echo "FASTQ_FILES: ${FASTQ_FILES[@]}"

# -------------------------
# Check array bounds
# -------------------------
NUM_FILES=${#FASTQ_FILES[@]}
if [ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]; then
    echo "Error: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID exceeds number of FASTQ files $NUM_FILES"
    exit 1
fi

# -------------------------
# Select file for this task
# -------------------------
FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
echo "Processing file: $FASTQ_FILE"

# -------------------------
# Determine mode
# -------------------------
if [[ $FASTQ_FILE == *"/AD_puro_only/"* ]]; then
    MODE="AD"
elif [[ $FASTQ_FILE == *"/RT_puro_only/"* ]]; then
    MODE="RT"
else
    echo "Error: Unable to determine mode for file $FASTQ_FILE"
    exit 1
fi
echo "Mode: $MODE"

# -------------------------
# Run TREBL Python script
# -------------------------
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/NKX2-2/scripts/trebl_experiment_pipeline_no_err_corr_puro_only_downsampling.py "$FASTQ_FILE" --mode "$MODE"