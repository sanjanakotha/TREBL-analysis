#!/bin/bash
#SBATCH --job-name=yeast_pool_C_umi_trebl_experiment_RT
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2
#SBATCH --array=0-20
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=24
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/logs/yeast_pool_C_umi_trebl_experiment_RT_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

echo "JOB STARTED on $(hostname) at $(date)"

# ---------- FILE LIST ----------
FILES=(
  /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/*fastq*
)

RT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Running on: $RT_FILE"

/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/yeast_pool_C_umi_trebl_experiment_RT_array.py "$RT_FILE"