Step 1 Protocol
===============

Overview
--------

This section describes how to use the TREBL analysis scripts to build and
refine an initial Step 1 map.

We will work through an example using the GCN4 Step 1 data. The data is
available on Savio, and you can start a new notebook without any special
conda environment. The initial map will likely take around 10 minutes
to create.

Steps
-----

1. Define barcode and domain structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, gather the preceding and post flanking sequences and the expected
length for each domain and activation domain (AD).

We will use :class:`~scripts.finder.Barcode` objects to organize barcode
names, preceder sequences, post sequences, and lengths.

.. code-block:: python

   import sys
   sys.path.append(
       "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/scripts"
   )
   import finder

   EC_AD = finder.Barcode(
       name="AD",
       preceder="GGCTAGC",
       post="",
       length=120,
   )

   EC_AD_BC = finder.Barcode(
       name="AD_BC",
       preceder="CGCGCC",
       post="",
       length=11,
   )

   EC_RPTR_BC = finder.Barcode(
       name="RPTR_BC",
       preceder="CTCGAG",
       post="",
       length=14,
   )

   bc_objects = [EC_AD, EC_AD_BC, EC_RPTR_BC]

2. Create a DuckDB database
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TREBL uses DuckDB in the background, which requires a database file.
Decide on a location where this file should be saved.

.. code-block:: python

   db_path = (
       "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/"
       "TREBL/duckdb/GCN4_final.db"
   )

3. Specify sequencing inputs and run the initial mapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will also need the following inputs:

- Sequencing data path(s) in ``.fastq`` or ``.fastq.gz`` format
- Design file path (if applicable)
- Whether to apply reverse complementation before searching with the
  provided flanking sequences

Pass these inputs to :class:`~scripts.initial_map.InitialMapper`.

.. code-block:: python

   step1_mapper = initial_map.InitialMapper(
       db_path=db_path,  # Where to save database output
       step_name="step1",  # Include time point / replicate
       seq_file=(
           "/global/scratch/projects/fc_mvslab/data/sequencing/"
           "CZB_Feb2024/A10_A11/results/"
           "A10_S1.fastq.gz.assembled.fastq"
       ),
       design_file_path=(
           "/global/scratch/projects/fc_mvslab/OpenProjects/EChase/"
           "TREBLEseq_ismaybethenewcibername/"
           "A10_sequencing/v2/current/a10_designfile.csv"
       ),
       bc_objects=bc_objects,
       reverse_complement=True,
   )

4. Refine map
~~~~~~~~~~~~~

After generating the initial Step 1 map, the next step is to refine it. 
We want to refine the map to only keep high-quality, high-confidence reads and to 
remove barcode collisions. 

The refinement process consists of several optional steps which can be performed 
in a flexible order. Each step is described below, with corresponding example
arguments which will be passed to :class:`~scripts.map_refiner.MapRefiner`.

Refinement steps
~~~~~~~~~~~~~~~~

Grouped
~~~~~~~

Groups reads and collapses duplicate reads to get counts per unique read. 
This is usually the first step in refinement.

Thresholded
~~~~~~~~~~~

Keeps reads above a certain abundance threshold. Must be run **after grouped**.

- If you know your read count threshold, set it directly:

  .. code-block:: python

     reads_threshold = 0

- If not, pass ``None``. You will be prompted to inspect the read histogram:

  .. code-block:: python

     reads_threshold = None

Barcode exists
~~~~~~~~~~~~~~

Remove rows with at least 1 empty barcode.

Quality
~~~~~~~

Filter reads based on sequence quality, meaning that domain and AD lengths must match expected length.

Unique target
~~~~~~~~~~~~~

- Which columns should be used for key–target matching? For example, if you want to check that one RPTR BC matches just one AD, you should pass the following argument:

  .. code-block:: python

     column_pairs = [("RPTR_BC", "AD")]

- But if you want to check that one RPTR BC matches a concatenated AD + AD BC, you should pass the following argument:

  .. code-block:: python

     column_pairs = [("RPTR_BC", ("AD", "AD_BC"))]

You can use any number of keys and/or targets, you just need to pass them in as a tuple.

- Do you want to allow flexibility in key/target matching? - For example, if the key is the RPTR barcode and the target is the AD, can a RPTR barcode map to multiple ADs as long as the most abundant AD makes up more than 90% of reads?

  .. code-block:: python

     min_fraction_major_target = 0.9  # RPTR can map to multiple ADs if >90% major

- Or should 1 RPTR barcode only map to 1 AD?

  .. code-block:: python

     min_fraction_major_target = 1

Designed
~~~~~~~~

Filters reads to keep only ADs that appear in the design file.

Recommended order
~~~~~~~~~~~~~~~~~

For consistency, we recommend applying steps in the following order:

1. grouped  
2. thresholded  
3. barcode_exists  
4. quality  
5. unique_target  
6. designed  

Thresholded **must come after grouped**, while the other steps can be applied 
flexibly.

We will pass the order explicitly:

.. code-block:: python

   map_order = [
       "grouped",
       "thresholded",
       "barcode_exists",
       "quality",
       "unique_target",
       "designed"
   ]

5. Refine the map
~~~~~~~~~~~~~~~~~

Finally, pass in desired arguments as above to :class:`~scripts.map_refiner.MapRefiner` constructor. 

Example: Apply all filters to keep sequences that:

- Have at least 50 reads
- Contain barcodes and ADs of correct length
- Map >90% of RPTR BC reads to the same AD
- Appear in the design file

.. code-block:: python

   refiner = map_refiner.MapRefiner(
       db_path=db_path,             # Same as before
       bc_objects=bc_objects,       # Same as before
       column_pairs=[("RPTR_BC", "AD")],
       reads_threshold=50,
       map_order=[
           "grouped",
           "thresholded",
           "quality",
           "unique_target",
           "designed",
       ],
       step_name="step1",
       output_figures_path="../../output/GCN4/figures/",  # Optional
   )

To run the refining:

.. code-block:: python

   refiner.refine_map_from_db()

Visualizing loss
----------------

Run this cell to create a bar graph showing total and unique reads at each refining step.

.. code-block:: python

   refiner.plot_loss()

Notes
-----

- In progress: Documentation for step 1 using UMI Tools to error correct reads.
