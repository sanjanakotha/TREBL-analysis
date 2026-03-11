from trebl_tools import pipelines, finder

pipeline = pipelines.TreblPipeline(
    db_path = "../../../duckdb/ChopTFs_pipeline.db",  
    design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/ChopTFDesign.csv",
    error_correction = False,   
    output_path = "../../../output/ChopTFs_pipeline"
)

AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_end = finder.Barcode(name="AD_end", preceder="", post="TGACTAG", length=40)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

step1_seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/GTAC_Jan2026/MZ/results/LIB129028_23FGWCLT3_S118_L006.fastq.gz.assembled.fastq" 

step1_map = pipeline.run_step_1(
    seq_file=step1_seq_file,         # FASTQ file input
    bc_objects=[AD, AD_end, AD_BC, RPTR_BC],           # Barcodes to map
    column_pairs=[("RPTR_BC", ("AD_end", "AD_BC"))],  # Check for collisions between reporter barcode and AD end + AD BC
    reads_threshold=50,              # Minimum number of reads to keep a barcode
    reverse_complement=True,        # Reverse complement reads for mapping
)