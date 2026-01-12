#!/bin/bash
#SBATCH --job-name=step1_err_corr
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/step1_err_corr_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

# Run the Python script
python $1