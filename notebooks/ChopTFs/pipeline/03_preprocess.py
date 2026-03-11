from trebl_tools import preprocess

# Filter AD reads
preprocess.run_fastp(input_dir="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Mar2026/TREBL_ChopTF_AD", 
                          output_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_AD_fastp")
preprocess.write_fastp_summary(output_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_AD_fastp")

# Filter RT reads
preprocess.run_fastp(input_dir="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Mar2026/TREBL_ChopTF_RP", 
                          output_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_RP_fastp")
preprocess.write_fastp_summary(output_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_RP_fastp")
