#!/bin/bash
#SBATCH --job-name=ChopTFs_step1
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/ChopTFs/ChopTFs_step1_%j.out

# Run the Python script
python $1