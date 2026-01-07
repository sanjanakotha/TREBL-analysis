#!/bin/bash
#SBATCH --job-name=umi_dedup
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/AD_1_15/umi_dedup_%j.out

# Run the Python script
python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/savio_jobs/umi_1_15_deduplication_test.py