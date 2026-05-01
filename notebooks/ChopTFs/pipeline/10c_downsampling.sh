#!/bin/bash
#SBATCH --job-name=trebl_downsample_ChopTFs
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:20:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/ChopTFs/pipeline/downsampling.sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --array=0,11,13,46,48,49

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
# FASTQ_FILES=(/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/NKX2-2_trebl_exp_chunks/AD_puro_only/*{50,100,200,400,600,800,1000}_chunks_part_{1..5}.fq.gz \
#              /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/NKX2-2_trebl_exp_chunks/RT_puro_only/*{50,100,200,400,600,800,1000}_chunks_part_{1..5}.fq.gz)

FASTQ_FILES=(/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_*_fastp/*_chunks/*{2,4,8,32,64,128,256,512}_chunks_part_{1..5}.fq.gz)

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
if [[ $FASTQ_FILE == *"/AD_chunks/"* ]]; then
    MODE="AD"
elif [[ $FASTQ_FILE == *"/RP_chunks/"* ]]; then
    MODE="RT"
else
    echo "Error: Unable to determine mode for file $FASTQ_FILE"
    exit 1
fi
echo "Mode: $MODE"

# -------------------------
# Run TREBL Python script
# -------------------------
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/ChopTFs/pipeline/10b_downsampling.py "$FASTQ_FILE" --mode "$MODE"