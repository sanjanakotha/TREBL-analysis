Examples
========

Here is an example workflow using TREBL-scripts.

Preliminary Mapping
-------------------
.. code-block:: python

    from mapping import BarcodeMapper

    # Intialize map object
    mapper = BarcodeMapper(seq_file = '/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A10_1_sequences.txt',
        design_file_path = "/global/scratch/projects/fc_mvslab/OpenProjects/EChase/TREBLEseq_ismaybethenewcibername/A10_sequencing/v2/current/a10_designfile.csv",
        bc_names = ["AD", "AD_BC", "RPTR_BC"],
        preceders = ["GGCTAGC", "CGCGCC", "CTCGAG"],
        posts = ["", "", ""],
        lengths = [120, 11, 14],
        reverse_complement = True)

    # Create map
    mapped_df = mapper.create_map()

    # Save map to parquet (slow for big datasets)
    mapped_df.save_parquet('/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/SK_CURRENT_A10_1_tbcRAW_v2.parquet')

Refining the Map
----------------
.. code-block:: python

    from map_refiner import MapRefiner

    # Initialize refiner object
    refiner = MapRefiner(db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/analysis_test.db",
        cols = ["AD", "AD_BC", "RPTR_BC"],
        reads_threshold = 50,
        column_pairs = [("AD", "AD_BC"), ("AD_BC", "RPTR_BC")]) # AD_BC maps to 1 AD, RPTR_BC matches to 1 AD_BC

    # Refine map starting with desired parquet file(s)
    refiner.refine_map_from_parquet('/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/SK_CURRENT_A10_1_tbcRAW_v2.parquet/*.parquet')

    # Save the final map
    refiner.save_map("map5", "final_map.csv")

    # Save a loss table
    loss_table = LossTable(db_path = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/duckdb/analysis_test.db",
                     cols = ["AD", "AD_BC", "RPTR_BC"])
