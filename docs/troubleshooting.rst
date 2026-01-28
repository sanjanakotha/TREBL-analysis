Troubleshooting
===============

This guide helps you resolve common issues when using TREBL Analysis.

Installation Issues
------------------

ModuleNotFoundError
~~~~~~~~~~~~~~~~~~~

**Problem**: ``ModuleNotFoundError: No module named 'scripts'`` or similar errors

**Solutions**:

1. **Add scripts to Python path**:

   .. code-block:: python

      import sys
      sys.path.append("/absolute/path/to/TREBL-analysis/scripts")

2. **Verify installation**:

   .. code-block:: bash

      pip list | grep -E "pandas|duckdb|dask"

3. **Reinstall dependencies**:

   .. code-block:: bash

      pip install -r docs/requirements.txt --force-reinstall

Python Version Errors
~~~~~~~~~~~~~~~~~~~~

**Problem**: ``SyntaxError`` or compatibility issues

**Solution**: TREBL requires Python 3.10+. Check your version:

.. code-block:: bash

   python --version

If needed, create a new environment:

.. code-block:: bash

   conda create -n trebl python=3.10
   conda activate trebl

Dependency Conflicts
~~~~~~~~~~~~~~~~~~~

**Problem**: Conflicting package versions

**Solution**: Use a fresh virtual environment:

.. code-block:: bash

   python -m venv trebl-env
   source trebl-env/bin/activate  # On Windows: trebl-env\Scripts\activate
   pip install -r docs/requirements.txt

Data Processing Issues
---------------------

No Barcodes Found
~~~~~~~~~~~~~~~~

**Problem**: Initial mapping finds 0 or very few barcodes

**Causes and Solutions**:

1. **Incorrect orientation**:
   
   Try toggling ``reverse_complement``:

   .. code-block:: python

      # If using True, try False
      mapper = initial_map.InitialMapper(
          # ... other parameters ...
          reverse_complement=False,  # Try opposite value
      )

2. **Wrong preceder sequences**:
   
   Verify your flanking sequences match the actual construct:

   .. code-block:: python

      # Double-check these match your plasmid design
      bc = finder.Barcode(
          name="AD",
          preceder="GGCTAGC",  # Verify this sequence
          post="",
          length=120,
      )

3. **Sequences not in FASTQ file**:
   
   Check a few reads manually:

   .. code-block:: bash

      zcat your_file.fastq.gz | head -20

4. **Case sensitivity**:
   
   Preceder sequences are case-sensitive. Ensure they match the case in your reads.

All Reads Filtered Out
~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Refinement removes all or most reads

**Solutions**:

1. **Lower read threshold**:

   .. code-block:: python

      refiner = map_refiner.MapRefiner(
          # ... other parameters ...
          reads_threshold=10,  # Try lower value
      )

2. **Remove strict filters**:
   
   Temporarily remove some filtering steps:

   .. code-block:: python

      map_order = [
          "grouped",
          "thresholded",
          # Comment out strict filters temporarily
          # "barcode_exists",
          # "quality",
      ]

3. **Check intermediate results**:
   
   Query the database after each step:

   .. code-block:: python

      import duckdb
      con = duckdb.connect(db_path)
      
      # List all tables
      tables = con.execute("SHOW TABLES").fetchall()
      print(tables)
      
      # Check counts at each step
      for table in tables:
          count = con.execute(f"SELECT COUNT(*) FROM {table[0]}").fetchone()
          print(f"{table[0]}: {count[0]} rows")

Barcode Length Mismatches
~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Quality filter removes reads due to length issues

**Solutions**:

1. **Check actual lengths**:

   .. code-block:: python

      import duckdb
      con = duckdb.connect(db_path)
      
      # Check length distribution
      lengths = con.execute("""
          SELECT LENGTH(AD) as ad_length, COUNT(*) as count
          FROM step1_initial
          GROUP BY ad_length
          ORDER BY count DESC
      """).fetchall()
      
      print(lengths)

2. **Adjust expected lengths**:
   
   Update barcode definitions to match actual data:

   .. code-block:: python

      AD = finder.Barcode(
          name="AD",
          preceder="GGCTAGC",
          post="",
          length=118,  # Adjust based on actual lengths
      )

Database Issues
--------------

Permission Denied
~~~~~~~~~~~~~~~~

**Problem**: Cannot create or write to database

**Solutions**:

1. **Check directory permissions**:

   .. code-block:: bash

      ls -ld /path/to/database/directory

2. **Use absolute path with write permissions**:

   .. code-block:: python

      import os
      db_path = os.path.abspath("./my_analysis.db")

3. **On HPC systems**, save to scratch or home directory:

   .. code-block:: python

      db_path = "/scratch/username/analysis.db"

Database Locked
~~~~~~~~~~~~~~

**Problem**: ``database is locked`` error

**Solutions**:

1. **Close other connections**:

   .. code-block:: python

      con.close()  # Close connection when done

2. **Check for other processes**:

   .. code-block:: bash

      # On Linux/macOS
      lsof my_analysis.db

3. **Delete lock file** (use with caution):

   .. code-block:: bash

      rm my_analysis.db.wal

Corrupted Database
~~~~~~~~~~~~~~~~~

**Problem**: Database appears corrupted or returns errors

**Solutions**:

1. **Start fresh**:

   .. code-block:: bash

      rm my_analysis.db
      # Re-run your analysis

2. **Check disk space**:

   .. code-block:: bash

      df -h

3. **Verify database integrity**:

   .. code-block:: python

      import duckdb
      con = duckdb.connect(db_path)
      con.execute("PRAGMA integrity_check").fetchall()

Memory Issues
------------

Out of Memory Errors
~~~~~~~~~~~~~~~~~~~

**Problem**: Process killed or ``MemoryError``

**Solutions**:

1. **Test with subset first**:

   .. code-block:: python

      mapper = initial_map.InitialMapper(
          # ... other parameters ...
          test_n_reads=10000,  # Process only 10k reads first
      )

2. **Process in batches**: Split your FASTQ file:

   .. code-block:: bash

      # Split into files with 1M reads each
      split -l 4000000 sequences.fastq sequences_part_

3. **Close unused connections**:

   .. code-block:: python

      con.close()

4. **Use a machine with more RAM**: Move to HPC cluster or larger instance

5. **Monitor memory usage**:

   .. code-block:: bash

      # On Linux
      htop
      # Or
      watch -n 1 free -h

File Format Issues
-----------------

FASTQ Format Errors
~~~~~~~~~~~~~~~~~~

**Problem**: Cannot read FASTQ file

**Solutions**:

1. **Verify FASTQ format**:

   .. code-block:: bash

      # Check first 20 lines
      head -20 sequences.fastq
      
      # For gzipped files
      zcat sequences.fastq.gz | head -20

   Each read should have 4 lines:
   
   - Line 1: ``@`` followed by sequence ID
   - Line 2: Sequence
   - Line 3: ``+``
   - Line 4: Quality scores

2. **Check file integrity**:

   .. code-block:: bash

      # For gzipped files
      gunzip -t sequences.fastq.gz

3. **Re-download or re-generate** if corrupted

Design File Issues
~~~~~~~~~~~~~~~~~

**Problem**: Design file not loaded correctly

**Solutions**:

1. **Verify CSV format**:
   
   Design file should be comma-separated with headers:

   .. code-block:: text

      AD,AD_BC,RPTR_BC
      ATCGATCG...,ATCGATCGATC,GCTAGCTAGCTAGC
      GCTAGCTA...,GCTAGCTAGCT,ATCGATCGATCGAT

2. **Check for special characters**:
   
   Ensure no extra commas, quotes, or special characters

3. **Use absolute path**:

   .. code-block:: python

      import os
      design_file = os.path.abspath("./design.csv")

Performance Issues
-----------------

Slow Processing
~~~~~~~~~~~~~~

**Problem**: Analysis takes too long

**Solutions**:

1. **Monitor progress**:
   
   Most functions include progress bars (via tqdm)

2. **Use SSD storage**: 
   
   Database I/O is faster on SSD vs HDD

3. **Check system resources**:

   .. code-block:: bash

      # CPU and memory usage
      top
      # Disk I/O
      iostat -x 1

4. **Optimize DuckDB settings**:

   .. code-block:: python

      con.execute("PRAGMA threads=8")  # Use multiple threads
      con.execute("PRAGMA memory_limit='16GB'")  # Increase memory

Visualization Issues
-------------------

Plots Not Displaying
~~~~~~~~~~~~~~~~~~~

**Problem**: ``refiner.plot_loss()`` doesn't show plots

**Solutions**:

1. **In Jupyter notebooks**:

   .. code-block:: python

      %matplotlib inline
      import matplotlib.pyplot as plt

2. **In scripts, save to file**:

   .. code-block:: python

      refiner.plot_loss()
      plt.savefig("loss_plot.png")
      plt.close()

3. **Use non-interactive backend**:

   .. code-block:: python

      import matplotlib
      matplotlib.use('Agg')  # Use before importing pyplot

Plotting Errors
~~~~~~~~~~~~~~

**Problem**: Matplotlib or seaborn errors

**Solutions**:

1. **Update packages**:

   .. code-block:: bash

      pip install --upgrade matplotlib seaborn

2. **Check backend**:

   .. code-block:: python

      import matplotlib
      print(matplotlib.get_backend())

HPC/Cluster Issues
-----------------

Job Failures on Slurm/SGE
~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Job killed or fails on cluster

**Solutions**:

1. **Request more memory**:

   .. code-block:: bash

      #SBATCH --mem=32G  # Increase memory allocation

2. **Request more time**:

   .. code-block:: bash

      #SBATCH --time=24:00:00  # Increase time limit

3. **Check job output**:

   .. code-block:: bash

      cat slurm-12345.out
      cat slurm-12345.err

4. **Test interactively first**:

   .. code-block:: bash

      srun --pty --mem=16G bash
      # Run your analysis interactively

Module Load Issues
~~~~~~~~~~~~~~~~~

**Problem**: Cannot load required modules

**Solutions**:

1. **List available modules**:

   .. code-block:: bash

      module avail python

2. **Load compatible Python**:

   .. code-block:: bash

      module load python/3.10

3. **Use conda instead**:

   .. code-block:: bash

      module load anaconda
      conda activate trebl

Error Correction Issues
----------------------

UMI Deduplication Errors
~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Error correction step fails

**Solutions**:

1. **Verify UMI format**: Ensure UMIs are correctly formatted

2. **Check error correction is enabled**:

   .. code-block:: python

      pipeline = TreblPipeline(
          # ... other parameters ...
          error_correction=True,  # Must be True
      )

3. **Ensure correct map order**:

   .. code-block:: python

      map_order = [
          "barcode_exists",
          "quality",
          "error_corrected",  # Error correction step
          "grouped",
          "thresholded",
          "unique_target",
          "designed"
      ]

Getting More Help
----------------

If your issue isn't covered here:

1. **Check existing documentation**:
   
   - :doc:`installation`
   - :doc:`quickstart`
   - :doc:`step1_protocol`
   - :doc:`scripts`

2. **Search GitHub Issues**:
   
   `<https://github.com/sanjanakotha/TREBL-analysis/issues>`_

3. **Open a new issue**:
   
   Include:
   
   - Python version (``python --version``)
   - Operating system
   - Complete error message
   - Minimal code to reproduce the issue
   - What you've already tried

4. **Check example notebooks**:
   
   The ``notebooks/`` directory contains working examples

Debugging Tips
-------------

Enable Verbose Output
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import logging
   logging.basicConfig(level=logging.DEBUG)

Check Database Contents
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import duckdb
   con = duckdb.connect(db_path)
   
   # List all tables
   print(con.execute("SHOW TABLES").fetchall())
   
   # Examine table structure
   print(con.execute("DESCRIBE step1_initial").fetchall())
   
   # Sample data
   print(con.execute("SELECT * FROM step1_initial LIMIT 5").fetchdf())

Test with Minimal Example
~~~~~~~~~~~~~~~~~~~~~~~~~

Create a small test FASTQ file:

.. code-block:: python

   test_fastq = """@read1
   GGCTAGCATCGATCGATCGCGCCCATCGATCGATCTCGAGATCGATCGATCGATCG
   +
   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
   @read2
   GGCTAGCATCGATCGATCGCGCCCATCGATCGATCTCGAGATCGATCGATCGATCG
   +
   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
   """
   
   with open("test.fastq", "w") as f:
       f.write(test_fastq)
   
   # Run analysis on test file
   # ...

Common Error Messages
--------------------

.. code-block:: text

   RuntimeError: VectorOperations::SimpleConstantCast

**Cause**: DuckDB type conversion error

**Solution**: Check that your design file columns match barcode names exactly

.. code-block:: text

   Catalog Error: Table with name "step1_initial" does not exist

**Cause**: Initial mapping hasn't been run or failed

**Solution**: Run ``mapper.run()`` before refinement

.. code-block:: text

   Binder Error: Referenced column "AD" not found

**Cause**: Barcode name doesn't match column name

**Solution**: Ensure barcode names in ``bc_objects`` match design file headers

Still Stuck?
-----------

Contact the developers:

- **GitHub**: `sanjanakotha/TREBL-analysis <https://github.com/sanjanakotha/TREBL-analysis>`_
- **Issues**: Open a detailed bug report with your error logs
