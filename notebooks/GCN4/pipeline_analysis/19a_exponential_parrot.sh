#!/bin/bash
# Job name:
#SBATCH --job-name=parrot_GCN4_exponential_more_seqs
#
# Account:
#SBATCH --account=fc_mvslab
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=4:00:00
#
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/pipeline_analysis/19a_exponential_parrot_chopTFs_merged.log
#
## Command(s) to run:
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot_fixed/bin/parrot-train \
#     '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/exponential_speed_input.txt' \
#     '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/exponential_speed' \
#     -d 'sequence' -c 1 --include-figs
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot_fixed/bin/parrot-train \
#     '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/exponential_strength_input_more_seqs.txt' \
#     '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/exponential_strength_more_seqs' \
#     -d 'sequence' -c 1 --include-figs
# /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot_fixed/bin/parrot-train \
#     '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/exponential_speed_input_k_more_seqs.txt' \
#     '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/exponential_speed_k_more_seqs' \
#     -d 'sequence' -c 1 --include-figs

/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot_fixed/bin/parrot-train \
    '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/chopTFs_gcn4_merged_input_strength.txt' \
    '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/chopTFs_gcn4_merged_strength' \
    -d 'sequence' -c 1 --include-figs
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot_fixed/bin/parrot-train \
    '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/chopTFs_gcn4_merged_input_k.txt' \
    '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/PARROT/chopTFs_gcn4_merged_k' \
    -d 'sequence' -c 1 --include-figs