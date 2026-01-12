#!/bin/bash
#SBATCH --job-name=fastp_array
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=32
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/notebooks/GCN4/last_fastp_file_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/umi_tools

fastp -i /global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2025/20250226_TREBL_MAZ06/MZ_EC_TREBL/RPTR_3_15_S18_R1_001.fastq.gz -o /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/RPTR_3_15_S18_R1_001_fastp.fastq.gz -w 32 --disable_adapter_trimming --json /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/logs/RPTR_3_15_S18_R1_001_fastp_report.json --html /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/logs/RPTR_3_15_S18_R1_001_fastp_report.html &> /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/logs/RPTR_3_15_S18_R1_001_fastp_report.out
