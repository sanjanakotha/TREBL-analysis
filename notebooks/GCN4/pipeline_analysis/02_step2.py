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

AD_step2_paths = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/LC_E1_step2/LC_E1_step2_testlibs_spikein/results/*AD*.assembled.fastq")
RPTR_step2_paths = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/LC_E1_step2/LC_E1_step2_testlibs_spikein/results/*RPTR*.assembled.fastq")

# Run Step 2 mapping
step2 = pipeline.run_step_2(
    AD_seq_file=AD_step2_paths,       # AD sequencing file
    AD_bc_objects=AD_bc_objects,         # AD barcodes
    RT_seq_file=RPTR_step2_paths,       # RT sequencing file
    RT_bc_objects=RT_bc_objects,         # RT barcodes
    reverse_complement=True,             # Search reverse complement
    reads_threshold_AD=10,               # Minimum reads for AD
    reads_threshold_RT=10,               # Minimum reads for RT
    step1_map_csv_path="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/step1.csv"  # Previous step 1 map in DuckDB
)
