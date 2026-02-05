#!/bin/bash
#SBATCH --job-name=trebl_experiment_pipeline_no_err_corr_rim
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/NKX2-2/logs/trebl_experiment_pipeline_no_err_corr_rim_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

echo "JOB STARTED on $(hostname) at $(date)"

# Run the Python script
python -u /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/NKX2-2/scripts/trebl_experiment_pipeline_no_err_corr_rim.py