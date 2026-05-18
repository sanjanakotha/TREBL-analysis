#!/bin/bash
#SBATCH --job-name=nardini_zscores
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3          
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          
#SBATCH --time=24:00:00
#SBATCH --output=21c_nardini.out

echo "Starting job on $(hostname)"

/global/home/users/sanjanakotha/.conda/envs/nardini_env/bin/python 21a_nardini_on_active.py

echo "Done"
