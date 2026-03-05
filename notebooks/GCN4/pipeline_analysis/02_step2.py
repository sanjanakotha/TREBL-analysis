from trebl_tools import pipelines, finder

pipeline = pipelines.TreblPipeline(
    db_path = "../../../duckdb/GCN4_final.db",  
    design_file_path = "../../../../../Marissa/DesignFiles/GCN4_Design.txt",
    error_correction = False,   
    output_path = "../../../output/GCN4_pipeline"
)
s
EC_AD = finder.Barcode(name = "AD",
                        preceder = "GGCTAGC",
                       post = "",
                       length = 120)

EC_AD_BC = finder.Barcode(name = "AD_BC",
                       preceder = "CGCGCC",
                       post = "",
                       length = 11)

EC_RPTR_BC = finder.Barcode(name = "RPTR_BC",
                       preceder = "CTCGAG",
                       post = "",
                       length = 14)

AD_step2_paths = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/LC_E1_step2/LC_E1_step2_testlibs_spikein/results/*AD*.assembled.fastq")
RPTR_step2_paths = glob.glob("/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/LC_E1_step2/LC_E1_step2_testlibs_spikein/results/*RPTR*.assembled.fastq")

# Produces histograms for AD and RT reads
# Helps pick reads_threshold_AD and reads_threshold_RT

# Run Step 2 mapping
step2 = pipeline.run_step_2(
    AD_seq_file=step2_AD_seq_file,       # AD sequencing file
    AD_bc_objects=AD_bc_objects,    s     # AD barcodes
    RT_seq_file=step2_RT_seq_file,       # RT sequencing file
    RT_bc_objects=RT_bc_objects,         # RT barcodes
    reverse_complement=True,             # Search reverse complement
    reads_threshold_AD=10,               # Minimum reads for AD
    reads_threshold_RT=10,               # Minimum reads for RT
    step1_map_csv_path="../../../output/GCN4_pipeline/step1.csv"  # Previous step 1 map in DuckDB
)

# Extract outputs
AD_step2 = step2["AD_step2"]                # Step 2 AD map
RT_step2 = step2["RT_step2"]                # Step 2 RT map
step1_step2_overlap = step2["step1_overlap"] # Overlap with Step 1 map