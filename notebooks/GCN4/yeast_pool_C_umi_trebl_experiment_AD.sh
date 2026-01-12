#!/bin/bash
#SBATCH --job-name=yeast_pool_C_umi_trebl_experiment_AD
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3_bigmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/yeast_pool_C_umi_trebl_experiment_AD_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

# Run the Python script
python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/yeast_pool_C_umi_trebl_experiment_AD.py