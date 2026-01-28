Using TREBL Pipelines
=====================

This guide focuses on using the high-level pipeline wrappers in ``pipelines.py``, which provide the simplest and most convenient way to run TREBL analyses.

Why Use Pipelines?
-----------------

The :class:`~scripts.pipelines.TreblPipeline` class wraps all the lower-level functions into an easy-to-use interface that:

- **Handles complexity**: Automatically manages database connections, intermediate tables, and data flow
- **Reduces errors**: Uses pre-configured, tested parameter combinations
- **Saves time**: One function call replaces many manual steps
- **Provides consistency**: Ensures analyses follow best practices
- **Simplifies batch processing**: Easily process multiple samples

When to Use Pipelines vs Lower-Level Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Use pipelines when:**

- You're running standard TREBL workflows
- You want quick, reliable results
- You're processing multiple samples
- You're new to TREBL analysis

**Use lower-level functions when:**

- You need fine-grained control over each step
- You're developing new analysis methods
- You have non-standard workflows

Getting Started with TreblPipeline
----------------------------------

Basic Setup
~~~~~~~~~~~

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   
   # Define your barcodes
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # Initialize the pipeline
   pipeline = TreblPipeline(
       db_path="./analysis.db",
       design_file_path="./design.csv",
       error_correction=False,  # Set True for UMI-based correction
       output_path="./output",   # Where to save results and figures
   )

Step 1 Analysis
--------------

Running Step 1 with Default Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to run a Step 1 analysis:

.. code-block:: python

   # Run Step 1 with default settings
   pipeline.run_step1(
       seq_file="./data/step1_sequences.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=50,
   )

This single function call:

1. Creates an initial map from your sequencing file
2. Applies the full refinement pipeline
3. Generates quality plots
4. Saves results to the output directory

Understanding run_step1 Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   pipeline.run_step1(
       seq_file,              # Path to FASTQ file
       bc_objects,            # List of Barcode objects
       reverse_complement,    # True/False - RC reads before mapping?
       reads_threshold,       # Minimum reads to keep a barcode combo
       column_pairs=None,     # Barcode-target relationships (optional)
       min_fraction_major_target=1,  # Collision threshold (optional)
   )

**seq_file** (str, required)
   Path to your FASTQ or FASTQ.gz sequencing file

**bc_objects** (list, required)
   List of :class:`~scripts.finder.Barcode` objects defining your construct

**reverse_complement** (bool, required)
   - ``True``: Apply reverse complement before barcode extraction
   - ``False``: Use sequences as-is
   - Most Illumina paired-end data needs ``True``

**reads_threshold** (int, required)
   Minimum number of reads required to keep a barcode combination
   - Lower values (10-20): More permissive, keeps more barcodes
   - Higher values (50-100): More stringent, only high-confidence barcodes

**column_pairs** (list, optional)
   Defines which barcodes should uniquely map to which targets. Default is ``[("RPTR_BC", "AD")]``
   
   Examples:
   
   .. code-block:: python
   
      # RPTR_BC should map to a single AD
      column_pairs=[("RPTR_BC", "AD")]
      
      # RPTR_BC should map to AD + AD_BC concatenation
      column_pairs=[("RPTR_BC", ("AD", "AD_BC"))]
      
      # Multiple relationships
      column_pairs=[
          ("RPTR_BC", "AD"),
          ("AD_BC", "AD")
      ]

**min_fraction_major_target** (float, optional)
   Controls how strict collision detection is. Default is 1 (no collisions allowed).
   
   - ``1.0``: Each key must map to exactly one target (strictest)
   - ``0.9``: Key can map to multiple targets if major one has ≥90% of reads
   - ``0.8``: Key can map to multiple targets if major one has ≥80% of reads

Complete Step 1 Example
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   
   # 1. Define barcodes
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # 2. Initialize pipeline
   pipeline = TreblPipeline(
       db_path="./gcn4_analysis.db",
       design_file_path="./gcn4_design.csv",
       error_correction=False,
       output_path="./gcn4_output",
   )
   
   # 3. Run Step 1
   pipeline.run_step1(
       seq_file="./data/gcn4_step1.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=50,
       column_pairs=[("RPTR_BC", "AD")],
       min_fraction_major_target=0.9,  # Allow 90% major
   )
   
   print("Step 1 complete! Check ./gcn4_output for results.")

Step 2 Analysis
--------------

Running Step 2
~~~~~~~~~~~~~

.. code-block:: python

   # Define Step 2 barcodes (may differ from Step 1)
   step2_bc_objects = [
       finder.Barcode(name="ADBC2", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="HawkBCs", preceder="TGGCAA", post="", length=20),
       finder.Barcode(name="RTBC", preceder="GGACGG", post="", length=14),
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
   ]
   
   # Run Step 2
   pipeline.run_step2(
       seq_file="./data/step2_sequences.fastq.gz",
       bc_objects=step2_bc_objects,
       reverse_complement=True,
       reads_threshold=20,
   )

Understanding run_step2 Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameters are similar to Step 1:

.. code-block:: python

   pipeline.run_step2(
       seq_file,              # Path to Step 2 FASTQ file
       bc_objects,            # List of Barcode objects (may differ from Step 1)
       reverse_complement,    # True/False
       reads_threshold,       # Minimum reads threshold
   )

Key differences from Step 1:

- Uses a different default refinement order (no ``unique_target`` step by default)
- Typically uses different barcode structures
- Often has lower read thresholds

Error Correction with Pipelines
-------------------------------

Enabling Error Correction
~~~~~~~~~~~~~~~~~~~~~~~~

To use UMI-based error correction:

.. code-block:: python

   # Initialize with error_correction=True
   pipeline = TreblPipeline(
       db_path="./analysis.db",
       design_file_path="./design.csv",
       error_correction=True,  # Enable error correction
       output_path="./output",
   )
   
   # Run normally - error correction happens automatically
   pipeline.run_step1(
       seq_file="./data/step1.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=50,
   )

What Error Correction Does
~~~~~~~~~~~~~~~~~~~~~~~~~~

When enabled, the pipeline:

1. Groups sequences by barcodes and UMIs
2. Corrects sequencing errors using a whitelist approach
3. Collapses corrected sequences
4. Filters based on corrected data

The error correction step is automatically inserted at the correct position in the refinement pipeline.

Batch Processing Multiple Samples
---------------------------------

Processing Multiple Step 1 Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the internal ``_run_initial_mappers`` method for batch processing:

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   
   # Define barcodes (same for all samples)
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # Initialize pipeline
   pipeline = TreblPipeline(
       db_path="./batch_analysis.db",
       design_file_path="./design.csv",
       output_path="./batch_output",
   )
   
   # Define all samples
   samples = [
       {
           "seq_file": "./data/sample1.fastq.gz",
           "step_name": "step1_sample1",
           "bc_objects": bc_objects,
           "reverse_complement": True,
           "design_file_path": "./design.csv",
       },
       {
           "seq_file": "./data/sample2.fastq.gz",
           "step_name": "step1_sample2",
           "bc_objects": bc_objects,
           "reverse_complement": True,
           "design_file_path": "./design.csv",
       },
       {
           "seq_file": "./data/sample3.fastq.gz",
           "step_name": "step1_sample3",
           "bc_objects": bc_objects,
           "reverse_complement": True,
           "design_file_path": "./design.csv",
       },
   ]
   
   # Run all initial mappings
   pipeline._run_initial_mappers(samples)
   
   # Then refine each one
   for i in range(1, 4):
       refiner = map_refiner.MapRefiner(
           db_path=pipeline.db_path,
           bc_objects=bc_objects,
           step_name=f"step1_sample{i}",
           reads_threshold=50,
           output_figures_path=str(pipeline.output_figures_path),
       )
       refiner.refine_map_from_db()

Advanced Pipeline Usage
----------------------

Testing with Subset of Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For quick testing or debugging:

.. code-block:: python

   # Initialize with test_n_reads
   pipeline = TreblPipeline(
       db_path="./test.db",
       design_file_path="./design.csv",
       test_n_reads=10000,  # Only process first 10,000 reads
       output_path="./test_output",
   )
   
   # Run normally - will automatically limit to 10k reads
   pipeline.run_step1(
       seq_file="./data/large_file.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=10,  # Use lower threshold for test data
   )

Custom Refinement Order
~~~~~~~~~~~~~~~~~~~~~~

By default, pipelines use pre-defined refinement orders. To customize:

.. code-block:: python

   from scripts import map_refiner
   
   # Run initial mapping using pipeline
   pipeline = TreblPipeline(
       db_path="./custom.db",
       design_file_path="./design.csv",
       output_path="./output",
   )
   
   # Use lower-level function for initial mapping
   mapper = initial_map.InitialMapper(
       db_path=pipeline.db_path,
       step_name="step1",
       seq_file="./data/sequences.fastq.gz",
       design_file_path=pipeline.design_file_path,
       bc_objects=bc_objects,
       reverse_complement=True,
   )
   mapper.run()
   
   # Custom refinement with your own order
   refiner = map_refiner.MapRefiner(
       db_path=pipeline.db_path,
       bc_objects=bc_objects,
       step_name="step1",
       reads_threshold=50,
       map_order=[
           "grouped",
           "thresholded",
           "quality",        # Skip some filters
           "designed",       # Or change order
       ],
   )
   refiner.refine_map_from_db()

Accessing Pipeline Results
--------------------------

Querying the Database
~~~~~~~~~~~~~~~~~~~~

After running the pipeline, access results via DuckDB:

.. code-block:: python

   import duckdb
   
   # Connect to pipeline database
   con = duckdb.connect(pipeline.db_path)
   
   # List all tables
   tables = con.execute("SHOW TABLES").fetchall()
   print("Available tables:", tables)
   
   # Get final refined data
   final_df = con.execute("""
       SELECT * FROM step1_AD_AD_BC_RPTR_BC_designed
   """).fetchdf()
   
   print(f"Final dataset: {len(final_df)} rows")
   
   # Summary statistics
   stats = con.execute("""
       SELECT 
           COUNT(*) as total_reads,
           COUNT(DISTINCT AD) as unique_ads,
           COUNT(DISTINCT RPTR_BC) as unique_rptr_bcs,
           AVG(reads) as avg_reads_per_combo
       FROM step1_AD_AD_BC_RPTR_BC_designed
   """).fetchone()
   
   print(f"Total reads: {stats[0]}")
   print(f"Unique ADs: {stats[1]}")
   print(f"Unique RPTR barcodes: {stats[2]}")
   print(f"Average reads per combination: {stats[3]:.1f}")
   
   con.close()

Exporting Results
~~~~~~~~~~~~~~~~

.. code-block:: python

   import duckdb
   import pandas as pd
   
   con = duckdb.connect(pipeline.db_path)
   
   # Export to CSV
   con.execute("""
       COPY step1_AD_AD_BC_RPTR_BC_designed 
       TO 'results.csv' (HEADER, DELIMITER ',')
   """)
   
   # Or export to pandas for further analysis
   df = con.execute("SELECT * FROM step1_AD_AD_BC_RPTR_BC_designed").fetchdf()
   
   # Analyze in pandas
   print(df.describe())
   
   # Export to Excel
   df.to_excel("results.xlsx", index=False)
   
   con.close()

Accessing Generated Figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline automatically saves figures to the output directory:

.. code-block:: python

   import os
   from pathlib import Path
   
   # List generated figures
   figures_dir = Path(pipeline.output_figures_path)
   figures = list(figures_dir.glob("*.png"))
   
   print("Generated figures:")
   for fig in figures:
       print(f"  - {fig.name}")

Real-World Examples
------------------

Example 1: GCN4 Analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   
   # GCN4 barcodes
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # Initialize
   pipeline = TreblPipeline(
       db_path="./gcn4.db",
       design_file_path="./gcn4_design.csv",
       output_path="./gcn4_results",
   )
   
   # Run Step 1
   pipeline.run_step1(
       seq_file="./data/A10_S1.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=50,
   )

Example 2: NKX2-2 with Error Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   
   # NKX2-2 barcodes
   bc_objects = [
       finder.Barcode(name="ADBC2", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="HawkBCs", preceder="TGGCAA", post="", length=20),
       finder.Barcode(name="RTBC", preceder="GGACGG", post="", length=14),
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
   ]
   
   # Initialize with error correction
   pipeline = TreblPipeline(
       db_path="./nkx22.db",
       design_file_path="./nkx22_design.csv",
       error_correction=True,  # Enable UMI correction
       output_path="./nkx22_results",
   )
   
   # Run Step 1
   pipeline.run_step1(
       seq_file="./data/nkx22_step1.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=20,
   )
   
   # Run Step 2
   pipeline.run_step2(
       seq_file="./data/nkx22_step2.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=15,
   )

Example 3: High-Throughput Batch Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   import glob
   
   # Define barcodes
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # Find all FASTQ files
   fastq_files = glob.glob("./data/batch_*.fastq.gz")
   
   # Process each file
   for i, fastq_file in enumerate(fastq_files, 1):
       print(f"Processing {fastq_file}...")
       
       # Create separate database for each sample
       db_path = f"./batch_sample{i}.db"
       
       pipeline = TreblPipeline(
           db_path=db_path,
           design_file_path="./design.csv",
           output_path=f"./batch_output/sample{i}",
       )
       
       pipeline.run_step1(
           seq_file=fastq_file,
           bc_objects=bc_objects,
           reverse_complement=True,
           reads_threshold=50,
       )
       
       print(f"Sample {i} complete!")

Pipeline Best Practices
----------------------

1. **Always test first**: Use ``test_n_reads`` with a small number before processing full datasets

2. **Use descriptive paths**: Name databases and output directories clearly:

   .. code-block:: python
   
      pipeline = TreblPipeline(
          db_path="./2025-01-gcn4-step1.db",  # Include date and experiment
          output_path="./results/gcn4/2025-01/",
      )

3. **Check output regularly**: Examine figures and summary statistics after each run

4. **Keep design files updated**: Ensure design files match your current library

5. **Document your parameters**: Save your pipeline script with comments:

   .. code-block:: python
   
      # GCN4 Step 1 Analysis - January 2025
      # Using 50-read threshold based on pilot data
      pipeline.run_step1(
          seq_file="./data/gcn4_jan2025.fastq.gz",
          reads_threshold=50,  # Based on pilot analysis
          # ... other parameters
      )

6. **Use version control**: Track your analysis scripts with git

7. **Separate raw and processed data**: Keep original FASTQ files separate from results

Common Pitfalls
--------------

**Forgetting to set reverse_complement correctly**
   Always check whether your sequencing reads need reverse complementing

**Using same database for different experiments**
   Use separate databases for each experiment to avoid conflicts

**Not checking intermediate results**
   Query the database after each major step to verify results

**Setting thresholds too high initially**
   Start with lower thresholds, examine results, then increase

**Ignoring generated figures**
   The loss plots show valuable information about filtering effectiveness

Next Steps
---------

- Learn about customization: :doc:`step1_protocol`
- Understand the underlying functions: :doc:`scripts`
- Troubleshoot issues: :doc:`troubleshooting`
- See working examples in the ``notebooks/`` directory

Summary
-------

The :class:`~scripts.pipelines.TreblPipeline` class provides:

- **Simple interface**: One function call for complete analyses
- **Sensible defaults**: Pre-configured parameters based on best practices
- **Automatic output**: Generates figures and tables automatically
- **Error correction**: Optional UMI-based correction
- **Batch processing**: Handle multiple samples efficiently

For most users, the pipeline wrappers provide everything needed for TREBL analysis!
