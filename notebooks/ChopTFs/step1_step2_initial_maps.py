import sys
sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/scripts")

import initial_map
import finder

# === Define Barcodes ===
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

# === Helper function to create, run, and preview a mapper ===
def run_mapper(db_path, seq_file, step_name, bc_objects, design_file_path=None, reverse_complement=True):
    mapper = initial_map.InitialMapper(
        db_path=db_path,
        seq_file=seq_file,
        step_name=step_name,
        bc_objects=bc_objects,
        reverse_complement=reverse_complement,
        design_file_path=design_file_path
    )
    mapper.create_map()
    return mapper

# === paths (same for all steps) ===
db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/ChopTFs.db"
design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/ChopTFDesign.csv"

# === Configurations for all steps ===
mapper_configs = [
    # STEP1
    dict(seq_file="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A11_S2.fastq.gz.assembled.fastq",
         step_name="step1",
         bc_objects=[AD, AD_BC, RPTR_BC],
         design_file_path=design_file_path),

    # STEP2 new
    dict(seq_file="/global/scratch/projects/fc_mvslab/data/sequencing/czb_new_sept2025/MAZ10/ChopTF/results/AD_Assembled/ChopTFstep2_2_AD_concat.fastq",
         step_name="step2_new",
         bc_objects=[AD, AD_BC],
         design_file_path=design_file_path),
    dict(seq_file="/global/scratch/projects/fc_mvslab/data/sequencing/czb_new_sept2025/MAZ10/ChopTF/results/RPTR_Assembled/ChopTFstep2_2_RPTR_concat.fastq",
         step_name="step2_new",
         bc_objects=[RPTR_BC],
         design_file_path=None),

    # STEP2 old
    dict(seq_file="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Jun2025/MAZ07/MZ_AU_SCU/results/AD_assembled/ChopTF_step2_AD_S1.fastq.gz.assembled.fastq",
         step_name="step2_old",
         bc_objects=[AD, AD_BC],
         design_file_path=design_file_path),
    dict(seq_file="/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Jun2025/MAZ07/MZ_AU_SCU/results/RPTR_assembled/ChopTF_step2_RPTR_S2.fastq.gz.assembled.fastq",
         step_name="step2_old",
         bc_objects=[RPTR_BC],
         design_file_path=None),
]

# === Run all mappers ===
mappers = [run_mapper(db_path=db_path, **cfg) for cfg in mapper_configs]
