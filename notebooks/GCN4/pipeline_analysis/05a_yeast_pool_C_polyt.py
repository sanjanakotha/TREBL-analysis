from trebl_tools import pipelines, finder
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

AD_bc_objects = [EC_AD, EC_AD_BC]
RT_bc_objects = [EC_RPTR_BC]

yeast_pool_C_PolyT_AD_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Sep2024/ciber2_iii_MZ001/EC_Ciber2_iii_Gcn4/results/assembled/AD/*")
yeast_pool_C_PolyT_RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Sep2024/ciber2_iii_MZ001/EC_Ciber2_iii_Gcn4/results/assembled/RPTR/*")

pipeline.trebl_experiment_analysis(
        AD_seq_files = yeast_pool_C_PolyT_AD_seq_files,
        AD_bc_objects = AD_bc_objects,
        RT_seq_files = yeast_pool_C_PolyT_RT_seq_files,
        RT_bc_objects = RT_bc_objects,
        reverse_complement = True,
        step1_map_csv_path="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/step1.csv",
        AD_umi_object=None,
        RT_umi_object=None,
        reads_threshold_AD=20,
        reads_threshold_RT=20,
        step_name_suffix="pool_C_PolyT_",
        umi_deduplication="simple")