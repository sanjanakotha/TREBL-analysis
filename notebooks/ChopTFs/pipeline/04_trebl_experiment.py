from trebl_tools import pipelines, finder, umi_deduplicate
import glob

pipeline = pipelines.TreblPipeline(
    db_path = "../../../duckdb/ChopTFs_pipeline.db",  
    design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/Gcn4_Design.csv",
    error_correction = False,   
    output_path = "../../../output/ChopTFs_pipeline"
)

# What sequences to search for?
AD_end = finder.Barcode(name="AD_end", preceder="", post="TGACTAG", length=40)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

AD_UMI = finder.Barcode(name = "UMI", preceder = "TGATTT", post = "", length = 12)
RT_UMI = finder.Barcode(name = "UMI", preceder = "TGTCAC", post = "", length = 12)

AD_bc_objects = [AD_end, AD_BC]
RT_bc_objects = [RPTR_BC]

AD_seq_files = []
#glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_AD_fastp/*")
RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/data/ChopTF/TREBL_ChopTF_RP_fastp/*")

pipeline.trebl_experiment_analysis(
        AD_seq_files = AD_seq_files,
        AD_bc_objects = AD_bc_objects,
        RT_seq_files = RT_seq_files,
        RT_bc_objects = RT_bc_objects,
        reverse_complement = True,
        step1_map_csv_path="../../../output/ChopTFs_pipeline/step1.csv",
        AD_umi_object=AD_UMI,
        RT_umi_object=RT_UMI,
        reads_threshold_AD=0,
        reads_threshold_RT=0,
        step_name_suffix="",
        umi_deduplication="both")