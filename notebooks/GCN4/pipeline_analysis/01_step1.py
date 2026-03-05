from trebl_tools import pipelines, finder

pipeline = pipelines.TreblPipeline(
    db_path = "../../../duckdb/GCN4_final.db",  
    design_file_path = "../../../../../Marissa/DesignFiles/GCN4_Design.txt",
    error_correction = False,   
    output_path = "../../../output/GCN4_pipeline"
)

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

step1_seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A10_S1.fastq.gz.assembled.fastq"


step1_map = pipeline.run_step_1(
    seq_file=step1_seq_file,         # FASTQ file input
    bc_objects=[EC_AD, EC_AD_BC, EC_RPTR_BC],           # Barcodes to map
    column_pairs=[("RPTR_BC", "AD")],  # Check for collisions between reporter barcode and AD
    reads_threshold=10,              # Minimum number of reads to keep a barcode
    reverse_complement=True         # Do not reverse complement reads for mapping
)