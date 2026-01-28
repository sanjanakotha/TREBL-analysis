Quick Start Guide
=================

This guide will help you run your first TREBL analysis in under 10 minutes using the **pipeline wrappers**.

.. note::
   
   **Recommended Approach**: Use the :class:`~scripts.pipelines.TreblPipeline` class for the simplest, most reliable workflow. For detailed pipeline documentation, see :doc:`pipelines_guide`.

What You'll Need
---------------

Before starting, ensure you have:

1. **TREBL Analysis installed**: See :doc:`installation`
2. **A FASTQ sequencing file**: Your raw sequencing data
3. **A design file** (CSV): Maps barcodes to their designed sequences
4. **Barcode information**: Flanking sequences and expected lengths

Recommended: Using the Pipeline Wrapper
---------------------------------------

The easiest way to run TREBL analysis is with the :class:`~scripts.pipelines.TreblPipeline` class.

Quick Example (3 Steps)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import sys
   sys.path.append("./scripts")
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
       db_path="./my_analysis.db",
       design_file_path="./design.csv",
       output_path="./output",  # Where to save results
   )
   
   # 3. Run complete Step 1 analysis
   pipeline.run_step1(
       seq_file="./data/sequences.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=50,
   )
   
   print("Analysis complete! Check ./output for results.")

That's it! The pipeline automatically:

- Extracts barcodes from sequencing reads
- Applies quality filters
- Removes barcode collisions
- Generates plots and statistics
- Saves all results

See :doc:`pipelines_guide` for comprehensive pipeline documentation.

Alternative: Step-by-Step Manual Approach
-----------------------------------------

For more control, you can use lower-level functions directly.

Step 0: Setup
~~~~~~~~~~~~~

.. code-block:: python

   import sys
   sys.path.append("./scripts")  # Adjust path as needed
   
   from scripts import finder
   from scripts import initial_map
   from scripts import map_refiner

Step 1: Define Your Barcodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create :class:`~scripts.finder.Barcode` objects for each domain in your construct:

.. code-block:: python

   # Define activation domain
   AD = finder.Barcode(
       name="AD",
       preceder="GGCTAGC",    # Sequence before the AD
       post="",               # Sequence after the AD
       length=120,            # Expected length in base pairs
   )
   
   # Define AD barcode
   AD_BC = finder.Barcode(
       name="AD_BC",
       preceder="CGCGCC",
       post="",
       length=11,
   )
   
   # Define reporter barcode
   RPTR_BC = finder.Barcode(
       name="RPTR_BC",
       preceder="CTCGAG",
       post="",
       length=14,
   )
   
   # Collect all barcodes
   bc_objects = [AD, AD_BC, RPTR_BC]

**Understanding Barcode Objects**:

- ``name``: Identifier for this barcode (used in output columns)
- ``preceder``: Sequence immediately before the barcode
- ``post``: Sequence immediately after the barcode (empty string if none)
- ``length``: Expected length of the barcode in base pairs

Step 2: Create Initial Map
~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the initial mapper to extract barcodes from your sequencing file:

.. code-block:: python

   # Specify file paths
   db_path = "./my_analysis.db"               # Where to save the database
   seq_file = "./data/sequences.fastq.gz"     # Your sequencing data
   design_file = "./data/design.csv"          # Design file (optional)
   
   # Create mapper
   mapper = initial_map.InitialMapper(
       db_path=db_path,
       step_name="step1",
       seq_file=seq_file,
       design_file_path=design_file,
       bc_objects=bc_objects,
       reverse_complement=True,  # Set to True if reads need RC
   )
   
   # Run mapping (this may take a few minutes)
   print("Starting initial mapping...")
   mapper.run()
   print("Initial mapping complete!")

**What happens during mapping?**

- Reads your FASTQ file
- Applies reverse complement if specified
- Searches for barcode preceder sequences
- Extracts barcodes between flanking sequences
- Stores results in a DuckDB database

Step 3: Refine the Map
~~~~~~~~~~~~~~~~~~~~~~

Apply quality filters to keep only high-confidence reads:

.. code-block:: python

   # Create refiner
   refiner = map_refiner.MapRefiner(
       db_path=db_path,
       bc_objects=bc_objects,
       step_name="step1",
       reads_threshold=50,                    # Minimum read count
       column_pairs=[("RPTR_BC", "AD")],     # Check RPTR maps to AD
       map_order=[
           "grouped",        # Group duplicate reads
           "thresholded",    # Apply read threshold
           "barcode_exists", # Remove empty barcodes
           "quality",        # Check length matches expected
           "unique_target",  # Remove barcode collisions
           "designed",       # Keep only designed sequences
       ],
   )
   
   # Run refinement
   print("Starting map refinement...")
   refiner.refine_map_from_db()
   print("Refinement complete!")

Step 4: Visualize Results
~~~~~~~~~~~~~~~~~~~~~~~~~

Plot the data loss at each filtering step:

.. code-block:: python

   # Generate loss plot
   refiner.plot_loss()

This creates a bar chart showing how many reads were filtered at each step.

Step 5: Access Results
~~~~~~~~~~~~~~~~~~~~~~

Retrieve your final filtered data:

.. code-block:: python

   # Get final table name
   final_table = refiner.get_final_table_name()
   print(f"Final table: {final_table}")
   
   # Query results
   import duckdb
   con = duckdb.connect(db_path)
   
   # Get summary statistics
   stats = con.execute(f"""
       SELECT COUNT(*) as total_reads,
              COUNT(DISTINCT AD) as unique_ads,
              COUNT(DISTINCT RPTR_BC) as unique_rptr_bcs
       FROM {final_table}
   """).fetchone()
   
   print(f"Total reads: {stats[0]}")
   print(f"Unique ADs: {stats[1]}")
   print(f"Unique RPTR barcodes: {stats[2]}")
   
   # Export to CSV
   df = con.execute(f"SELECT * FROM {final_table}").fetchdf()
   df.to_csv("results.csv", index=False)
   
   con.close()

Complete Example Script (Manual Approach)
-----------------------------------------

Here's a complete script using lower-level functions:

.. code-block:: python

   import sys
   sys.path.append("./scripts")
   from scripts import finder, initial_map, map_refiner
   import duckdb
   
   # 1. Define barcodes
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # 2. File paths
   db_path = "./analysis.db"
   seq_file = "./data/sequences.fastq.gz"
   design_file = "./data/design.csv"
   
   # 3. Initial mapping
   mapper = initial_map.InitialMapper(
       db_path=db_path,
       step_name="step1",
       seq_file=seq_file,
       design_file_path=design_file,
       bc_objects=bc_objects,
       reverse_complement=True,
   )
   mapper.run()
   
   # 4. Refinement
   refiner = map_refiner.MapRefiner(
       db_path=db_path,
       bc_objects=bc_objects,
       step_name="step1",
       reads_threshold=50,
       column_pairs=[("RPTR_BC", "AD")],
       map_order=[
           "grouped", "thresholded", "barcode_exists",
           "quality", "unique_target", "designed"
       ],
   )
   refiner.refine_map_from_db()
   
   # 5. Results
   refiner.plot_loss()
   
   # Export results
   con = duckdb.connect(db_path)
   final_table = refiner.get_final_table_name()
   df = con.execute(f"SELECT * FROM {final_table}").fetchdf()
   df.to_csv("results.csv", index=False)
   con.close()
   
   print("Analysis complete! Results saved to results.csv")

Which Approach Should I Use?
----------------------------

**Use the Pipeline Wrapper** (:doc:`pipelines_guide`) when:

- You want the simplest, fastest workflow
- You're running standard TREBL analyses
- You're new to TREBL analysis
- You want automatic figure generation

**Use Manual Steps** (shown above) when:

- You need fine-grained control over each step
- You're developing new analysis methods
- You have non-standard workflows

.. tip::
   
   **For most users, we recommend starting with the pipeline wrapper!** See :doc:`pipelines_guide` for complete documentation.

Analyzing Multiple Samples
--------------------------

Using the pipeline for batch processing:

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   import glob
   
   # Find all FASTQ files
   fastq_files = glob.glob("./data/*.fastq.gz")
   
   # Process each file
   for fastq_file in fastq_files:
       pipeline = TreblPipeline(
           db_path=f"./{fastq_file.stem}.db",
           design_file_path="./design.csv",
           output_path=f"./output/{fastq_file.stem}",
       )
       
       pipeline.run_step1(
           seq_file=fastq_file,
           bc_objects=bc_objects,
           reverse_complement=True,
           reads_threshold=50,
       )

See :doc:`pipelines_guide` for more batch processing examples.

Next Steps
---------

Now that you've completed your first analysis:

1. **Understand the protocol**: Read :doc:`step1_protocol` for detailed explanations
2. **Explore refinement options**: Learn about different filtering strategies
3. **Review the API**: Check :doc:`scripts` for all available functions
4. **Analyze your own data**: Adapt these examples to your experiments

Common Parameters
----------------

Understanding Key Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**reverse_complement**
   Set to ``True`` if your sequencing reads are in the opposite orientation from your barcode definitions. Most Illumina paired-end data requires this.

**reads_threshold**
   Minimum number of reads required to keep a barcode combination. Higher values = more stringent filtering. Typical values: 10-100.

**column_pairs**
   Specifies which barcodes should uniquely map to which targets. For example:
   
   - ``[("RPTR_BC", "AD")]``: Each reporter barcode should map to one AD
   - ``[("RPTR_BC", ("AD", "AD_BC"))]``: RPTR maps to AD+AD_BC concatenation

**map_order**
   Order of filtering steps. Must include ``"grouped"`` before ``"thresholded"``. Other steps are flexible.

Typical Refinement Orders
~~~~~~~~~~~~~~~~~~~~~~~~~

**Standard workflow** (no error correction):

.. code-block:: python

   map_order = [
       "grouped",
       "thresholded",
       "barcode_exists",
       "quality",
       "unique_target",
       "designed"
   ]

**With error correction**:

.. code-block:: python

   map_order = [
       "barcode_exists",
       "quality",
       "error_corrected",
       "grouped",
       "thresholded",
       "unique_target",
       "designed"
   ]

Tips and Best Practices
-----------------------

1. **Start with test data**: Use ``test_n_reads=1000`` to quickly test your pipeline before running on full datasets

2. **Check intermediate results**: Query the database between steps to verify results

3. **Adjust thresholds iteratively**: Start with low thresholds, examine results, then increase

4. **Save your scripts**: Keep a record of the exact parameters you used

5. **Use absolute paths**: Especially on HPC systems, use absolute paths for files and databases

6. **Monitor memory**: Large FASTQ files can consume significant memory

Troubleshooting
--------------

**No barcodes found**
   - Check that ``preceder`` sequences are correct
   - Try both ``reverse_complement=True`` and ``reverse_complement=False``
   - Verify your FASTQ file is formatted correctly

**Too few reads after filtering**
   - Lower your ``reads_threshold``
   - Remove some filtering steps from ``map_order``
   - Check quality of input data

**Database errors**
   - Ensure directory exists and you have write permissions
   - Use absolute paths instead of relative paths
   - Check disk space

**Memory errors**
   - Process fewer reads at once using ``test_n_reads``
   - Close other applications
   - Use a machine with more RAM

Getting Help
-----------

- **Documentation**: :doc:`index`
- **GitHub Issues**: `Report bugs or ask questions <https://github.com/sanjanakotha/TREBL-analysis/issues>`_
- **Examples**: Check the ``notebooks/`` directory for Jupyter notebook examples
