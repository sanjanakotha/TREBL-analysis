from trebl_tools import pipelines, finder, umi_deduplicate
import glob

pipeline = pipelines.TreblPipeline(
    db_path = "../../../duckdb/GCN4_pipeline.db",  
    design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/Gcn4_Design.csv",
    error_correction = False,   
    output_path = "../../../output/GCN4_pipeline"
)

EC_AD = finder.Barcode(name = "AD", preceder = "GGCTAGC", post = "", length = 120)
EC_AD_BC = finder.Barcode(name = "AD_BC", preceder = "CGCGCC", post = "", length = 11)
EC_RPTR_BC = finder.Barcode(name = "RPTR_BC", preceder = "CTCGAG", post = "", length = 14)

AD_UMI = finder.Barcode(name = "UMI", preceder = "TGATTT", post = "", length = 12)
RT_UMI = finder.Barcode(name = "UMI", preceder = "TGTCAC", post = "", length = 12)

AD_bc_objects = [EC_AD, EC_AD_BC]
RT_bc_objects = [EC_RPTR_BC]

## Running fastp on single end RT reads
# preprocess.run_fastp(input_dir="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2025/20250226_TREBL_MAZ06/MZ_EC_TREBL/", 
#                           output_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/", 
#                           script_path="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/savio_jobs/fastp.sh")

yeast_pool_C_umi_AD_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/20250218_MZCCSCU_MedGenome/MZ/results/assembled/*")
yeast_pool_C_umi_RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/GCN4_pool_C_UMI_RPTR_fastp/*")

pipeline.trebl_experiment_analysis(
        AD_seq_files = yeast_pool_C_umi_AD_seq_files,
        AD_bc_objects = AD_bc_objects,
        RT_seq_files = yeast_pool_C_umi_RT_seq_files,
        RT_bc_objects = RT_bc_objects,
        reverse_complement = True,
        step1_map_csv_path="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/step1.csv",
        AD_umi_object=AD_UMI,
        RT_umi_object=RT_UMI,
        reads_threshold_AD=0,
        reads_threshold_RT=0,
        step_name_suffix="pool_C_umi_",
        umi_deduplication="both")