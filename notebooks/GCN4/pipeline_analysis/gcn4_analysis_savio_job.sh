#!/bin/bash
#SBATCH --job-name=gcn4_analysis_redone
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio2_bigmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/pipeline_analysis.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

# echo "Starting step1 analysis"
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/01_step1.py

# echo "Starting step2 analysis"
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/02_step2.py

# echo "Starting yeast pool A analysis"
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/03_yeast_pool_A.py

# echo "Starting yeast pool B analysis"
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/04_yeast_pool_B.py

# echo "Start yeast pool C PolyT analysis"
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/05_yeast_pool_C_polyt.py

echo "Start yeast pool C UMI analysis"
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/06_yeast_pool_C_umi.py