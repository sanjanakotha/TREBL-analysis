#!/bin/bash
#SBATCH --job-name=GCN4_step1_err_corr
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/GCN4_step1_err_corr_%j.out

# Run the Python script
python $1