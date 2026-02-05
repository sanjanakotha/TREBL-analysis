#!/bin/bash
#SBATCH --job-name=yeast_pool_C_umi_trebl_experiment_RT
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/logs/yeast_pool_C_umi_trebl_experiment_RT_%j.out
#SBATCH --error=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/logs/yeast_pool_C_umi_trebl_experiment_RT_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

echo "JOB STARTED on $(hostname) at $(date)"

source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env

# Run the Python script
python -u /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/yeast_pool_C_umi_trebl_experiment_RT.py