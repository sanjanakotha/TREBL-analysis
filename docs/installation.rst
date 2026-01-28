Installation Guide
==================

This guide will help you install TREBL Analysis and its dependencies.

Prerequisites
-------------

Before installing TREBL Analysis, ensure you have:

- **Python 3.10 or higher**: Check your version with ``python --version`` or ``python3 --version``
- **pip**: Python package installer (usually comes with Python)
- **Git**: For cloning the repository

Recommended System Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Memory**: At least 8GB RAM (16GB+ recommended for large datasets)
- **Storage**: Sufficient disk space for sequencing files and DuckDB databases
- **OS**: Linux, macOS, or Windows (Linux/macOS recommended)

Installation Methods
-------------------

Method 1: Standard Installation (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Clone the repository**:

   .. code-block:: bash

      git clone https://github.com/sanjanakotha/TREBL-analysis.git
      cd TREBL-analysis

2. **Install dependencies**:

   .. code-block:: bash

      pip install -r docs/requirements.txt

3. **Verify installation**:

   .. code-block:: python

      import sys
      sys.path.append("./scripts")
      from scripts import finder, initial_map, map_refiner
      print("TREBL Analysis successfully installed!")

Method 2: Virtual Environment (Recommended for Isolation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using a virtual environment keeps TREBL dependencies separate from your system Python.

**Using venv**:

.. code-block:: bash

   # Clone repository
   git clone https://github.com/sanjanakotha/TREBL-analysis.git
   cd TREBL-analysis

   # Create virtual environment
   python3 -m venv trebl-env

   # Activate virtual environment
   # On Linux/macOS:
   source trebl-env/bin/activate
   # On Windows:
   trebl-env\Scripts\activate

   # Install dependencies
   pip install -r docs/requirements.txt

**Using conda**:

.. code-block:: bash

   # Clone repository
   git clone https://github.com/sanjanakotha/TREBL-analysis.git
   cd TREBL-analysis

   # Create conda environment
   conda create -n trebl python=3.10
   conda activate trebl

   # Install dependencies
   pip install -r docs/requirements.txt

Method 3: Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you plan to modify the code or contribute:

.. code-block:: bash

   git clone https://github.com/YOUR_USERNAME/TREBL-analysis.git
   cd TREBL-analysis
   
   # Add upstream remote
   git remote add upstream https://github.com/sanjanakotha/TREBL-analysis.git
   
   # Install dependencies
   pip install -r docs/requirements.txt

Dependencies
-----------

TREBL Analysis requires the following Python packages:

Core Dependencies
~~~~~~~~~~~~~~~~~

- **pandas**: Data manipulation and analysis
- **dask**: Parallel computing with task scheduling
- **duckdb**: Embedded SQL database for efficient data processing
- **biopython**: Biological computation tools
- **numpy**: Numerical computing

Visualization
~~~~~~~~~~~~~

- **matplotlib**: Plotting library
- **seaborn**: Statistical data visualization

Utilities
~~~~~~~~~

- **pyarrow**: Columnar data format (required by DuckDB)
- **tqdm**: Progress bars

Documentation (Optional)
~~~~~~~~~~~~~~~~~~~~~~~

Only needed if building documentation locally:

- **sphinx**: Documentation generator
- **furo**: Sphinx theme

HPC/Cluster Installation
-----------------------

If you're working on an HPC cluster (e.g., Savio):

1. **Load required modules** (if available):

   .. code-block:: bash

      module load python/3.10
      module load git

2. **Create isolated environment**:

   .. code-block:: bash

      # Clone to your scratch or project directory
      cd $SCRATCH  # or your preferred location
      git clone https://github.com/sanjanakotha/TREBL-analysis.git
      cd TREBL-analysis

3. **Install to user directory**:

   .. code-block:: bash

      pip install --user -r docs/requirements.txt

4. **Add to Python path** in your job scripts:

   .. code-block:: python

      import sys
      sys.path.append("/path/to/TREBL-analysis/scripts")

Verifying Your Installation
---------------------------

Test Basic Import
~~~~~~~~~~~~~~~~~

.. code-block:: python

   import sys
   sys.path.append("./scripts")
   
   from scripts import finder
   from scripts import initial_map
   from scripts import map_refiner
   from scripts import pipelines
   
   print("All modules imported successfully!")

Test DuckDB Connection
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import duckdb
   
   # Create test database
   con = duckdb.connect("test.db")
   con.execute("CREATE TABLE test (id INTEGER, name VARCHAR)")
   con.execute("INSERT INTO test VALUES (1, 'hello')")
   result = con.execute("SELECT * FROM test").fetchall()
   print(f"DuckDB test: {result}")
   con.close()
   
   # Clean up
   import os
   os.remove("test.db")

Test Barcode Creation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import sys
   sys.path.append("./scripts")
   from scripts import finder
   
   # Create a test barcode
   test_bc = finder.Barcode(
       name="TEST",
       preceder="ATCG",
       post="GCTA",
       length=10
   )
   
   print(f"Barcode created: {test_bc.name}")
   print(f"Preceder: {test_bc.preceder}")
   print(f"Length: {test_bc.length}")

Troubleshooting
--------------

Common Issues
~~~~~~~~~~~~~

**Import Error: Module not found**

If you get ``ModuleNotFoundError``, ensure:

- You've installed all dependencies: ``pip install -r docs/requirements.txt``
- You've added the scripts directory to your Python path:

  .. code-block:: python

     import sys
     sys.path.append("/absolute/path/to/TREBL-analysis/scripts")

**DuckDB Connection Error**

If you can't connect to DuckDB:

- Ensure you have write permissions in the directory where you're creating the database
- Check that the path is valid and the directory exists
- Try using an absolute path instead of a relative path

**Memory Error with Large Files**

For large sequencing files:

- Ensure you have sufficient RAM (16GB+ recommended)
- Consider using the ``test_n_reads`` parameter to process a subset first
- Close other memory-intensive applications

**Python Version Incompatibility**

TREBL requires Python 3.10+. If you have an older version:

- Install Python 3.10 or higher
- Use conda to create an environment with the correct version:

  .. code-block:: bash

     conda create -n trebl python=3.10
     conda activate trebl

**Permission Errors on HPC**

On shared computing systems:

- Use ``pip install --user`` to install to your user directory
- Or create a virtual environment in your home/scratch directory

Updating TREBL Analysis
-----------------------

To update to the latest version:

.. code-block:: bash

   cd TREBL-analysis
   git pull origin main
   pip install -r docs/requirements.txt --upgrade

Uninstalling
-----------

To remove TREBL Analysis:

1. **Delete the repository**:

   .. code-block:: bash

      rm -rf /path/to/TREBL-analysis

2. **Remove dependencies** (if using virtual environment):

   .. code-block:: bash

      # For venv
      deactivate
      rm -rf trebl-env

      # For conda
      conda deactivate
      conda env remove -n trebl

Getting Help
-----------

If you encounter issues:

1. Check this troubleshooting section
2. Review the `Quick Start Guide <quickstart.html>`_
3. Search existing `GitHub Issues <https://github.com/sanjanakotha/TREBL-analysis/issues>`_
4. Open a new issue with:
   - Your Python version
   - Operating system
   - Full error message
   - Steps to reproduce

Next Steps
----------

After installation, proceed to:

- :doc:`quickstart` - Get started with your first analysis
- :doc:`step1_protocol` - Detailed Step 1 protocol
- :doc:`scripts` - API reference documentation
