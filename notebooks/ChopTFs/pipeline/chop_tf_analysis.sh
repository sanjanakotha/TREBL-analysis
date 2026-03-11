#!/bin/bash
#SBATCH --job-name=choptf_analysis
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/ChopTFs/pipeline/pipeline_analysis.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

# echo "Starting step1 analysis"
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/ChopTFs/pipeline/01_step1.py

echo "Starting trebl experiment analysis"
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/ChopTFs/pipeline/04_trebl_experiment.py