TREBL Analysis Documentation
============================

Welcome to the TREBL Analysis documentation!

TREBL (Transcription factor REporter Barcode Library) is a high-throughput sequencing-based method for mapping transcription factor activation domains to DNA barcodes. This documentation provides comprehensive guides, tutorials, and API references for analyzing TREBL sequencing data.

.. note::
   
   **New to TREBL Analysis?** 
   
   1. Start with :doc:`installation` to set up the tools
   2. Follow :doc:`quickstart` for your first analysis
   3. Learn the recommended workflow in :doc:`pipelines_guide`

Getting Started
--------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart
   pipelines_guide

Protocols & Guides
-----------------

Detailed workflows and protocols for TREBL experiments:

.. toctree::
   :maxdepth: 2
   :caption: Protocols

   step1_protocol

.. note::
   
   Step 2 and TREBL experiment protocols are coming soon!

API Reference
-------------

Complete documentation of all modules and functions:

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   scripts

Help & Support
-------------

.. toctree::
   :maxdepth: 2
   :caption: Help

   troubleshooting

Quick Links
----------

* **Installation**: :doc:`installation`
* **Quick Start**: :doc:`quickstart`
* **Pipeline Guide** (Recommended): :doc:`pipelines_guide`
* **Step 1 Protocol**: :doc:`step1_protocol`
* **Troubleshooting**: :doc:`troubleshooting`
* **GitHub Repository**: `sanjanakotha/TREBL-analysis <https://github.com/sanjanakotha/TREBL-analysis>`_

Recommended Workflow
-------------------

For most users, we recommend using the :class:`~scripts.pipelines.TreblPipeline` wrapper class:

.. code-block:: python

   from scripts.pipelines import TreblPipeline
   from scripts import finder
   
   # Define barcodes
   bc_objects = [
       finder.Barcode(name="AD", preceder="GGCTAGC", post="", length=120),
       finder.Barcode(name="AD_BC", preceder="CGCGCC", post="", length=11),
       finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="", length=14),
   ]
   
   # Initialize and run
   pipeline = TreblPipeline(
       db_path="./analysis.db",
       design_file_path="./design.csv",
       output_path="./output",
   )
   
   pipeline.run_step1(
       seq_file="./data/sequences.fastq.gz",
       bc_objects=bc_objects,
       reverse_complement=True,
       reads_threshold=50,
   )

See :doc:`pipelines_guide` for detailed documentation on using the pipeline wrappers.

Workflow Overview
----------------

A typical TREBL analysis consists of three main stages:

**Step 1: Initial Mapping**
   Extract barcodes from sequencing reads and create an initial map linking barcodes to activation domains.

**Step 2: Validation** (coming soon)
   Validate the barcode-to-AD mappings through additional sequencing.

**Step 3: TREBL Experiment** (coming soon)
   Perform the full TREBL experiment and analyze results.

Key Features
-----------

* **Simple pipeline interface**: One function call for complete analyses
* **Flexible barcode extraction**: Define custom barcode structures with flanking sequences
* **Efficient processing**: Uses DuckDB for fast analysis of large datasets
* **Quality control**: Multi-step filtering pipeline to ensure data quality
* **Error correction**: Optional UMI-based error correction
* **Visualization**: Built-in plotting functions for quality metrics
* **Batch processing**: Handle multiple samples in a single run

Indices and Tables
-----------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
