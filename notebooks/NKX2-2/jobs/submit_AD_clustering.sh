#!/bin/bash
#SBATCH --job-name=NKX2-2_AD_clustering
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/NKX2-2/NKX2-2_AD_clustering_%j.out

source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/umi_tools

python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/NKX2-2/UMI_Clusterer_NKX2-2.py