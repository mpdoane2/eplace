Installation
============

This guide covers the installation of ePLACE and all its dependencies.

Prerequisites
-------------

Python Version
~~~~~~~~~~~~~~

ePLACE requires Python 3.8 or higher.

.. code-block:: bash

   # Check your Python version
   python --version

Using Conda/Mamba (Recommended)
--------------------------------

1. Create a conda environment for ePLACE:

.. code-block:: bash

   mamba create -yn eplace 'python>=3.12'
   conda activate eplace

2. Install the required dependencies:

.. code-block:: bash

   mamba install -y bioconda::blast bioconda::pytaxonkit bioconda::iqtree bioconda::mafft

.. note::
   At the time of writing there is an `issue <https://github.com/bioforensics/pytaxonkit/issues/50>`_ with conda not installing the
   most current version of pytaxonkit if you are using python >=3.12. This code works with older versions of pytaxonkit.

.. note::
   You will need to download and set up the NCBI taxonomy databases for pytaxonkit; see the 
   `taxonkit documentation <https://bioinf.shenwei.me/taxonkit/usage/#before-use>`_ for detailed 
   instructions of which NCBI taxonomy files to download.

3. Install ePLACE:

.. code-block:: bash

   pip install git+https://github.com/linsalrob/eplace.git

Standard Installation
---------------------

From GitHub
~~~~~~~~~~~

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/linsalrob/eplace.git
   cd eplace

   # Install the package
   pip install .

   # Or install in development mode
   pip install -e .

   # Or with development dependencies
   pip install -e ".[dev]"

After installation, the ``eplace`` command will be available in your environment.

External Dependencies
---------------------

ePLACE requires several external bioinformatics tools to be installed separately.

BLAST+ Tools (Required)
~~~~~~~~~~~~~~~~~~~~~~~~

BLAST+ tools (blastn, blastdbcmd) must be installed:

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get install ncbi-blast+

**macOS with Homebrew:**

.. code-block:: bash

   brew install blast

**From Source:**

Download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

TaxonKit (Required)
~~~~~~~~~~~~~~~~~~~

TaxonKit is required for taxonomy lookup:

**Via Conda:**

.. code-block:: bash

   conda install -c bioconda taxonkit

**From Binary:**

Download from: https://github.com/shenwei356/taxonkit/releases

MAFFT (Optional)
~~~~~~~~~~~~~~~~

MAFFT is required for sequence alignment (can be skipped with ``--skip-alignment`` flag):

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get install mafft

**macOS with Homebrew:**

.. code-block:: bash

   brew install mafft

**Via Conda:**

.. code-block:: bash

   conda install -c bioconda mafft

IQTree (Optional)
~~~~~~~~~~~~~~~~~

IQTree is required for phylogenetic tree building (can be skipped with ``--skip-alignment`` flag):

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get install iqtree

**macOS with Homebrew:**

.. code-block:: bash

   brew install iqtree

**Via Conda:**

.. code-block:: bash

   conda install -c bioconda iqtree

Verify Installation
-------------------

After installation, verify that ePLACE is installed correctly:

.. code-block:: bash

   # Check the installation
   eplace --help

   # Show version
   eplace --version

   # Test that BLAST is available
   blastn -version

   # Test that TaxonKit is available
   taxonkit version

Download NCBI Database
----------------------

Before running analyses, you need to download the NCBI BLAST database:

.. code-block:: bash

   # Download core_nt database to default location
   eplace download

   # Force redownload even if database exists
   eplace download --force

The database will be stored in ``$BLASTDB`` if set, otherwise in ``~/blastdb``.

.. warning::
   The NCBI core_nt database is very large (several GB). Ensure you have sufficient disk space 
   and bandwidth before downloading. The download may take a significant amount of time depending 
   on your internet connection.

Development Installation
------------------------

For developers who want to contribute to ePLACE:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/linsalrob/eplace.git
   cd eplace

   # Install in development mode with dev dependencies
   pip install -e ".[dev]"

   # Run tests
   pytest tests/ -v

   # Run tests with coverage
   pytest tests/ --cov=eplace_lib --cov-report=html

Troubleshooting
---------------

Command not found
~~~~~~~~~~~~~~~~~

If ``eplace`` command is not found after installation:

.. code-block:: bash

   # Check if it's installed
   pip show eplace

   # Reinstall
   pip install --force-reinstall .

   # Or add to PATH
   export PATH="$HOME/.local/bin:$PATH"

Permission denied
~~~~~~~~~~~~~~~~~

If you encounter permission errors during installation:

.. code-block:: bash

   # Install for current user only
   pip install --user .

BLASTDB environment variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To set a custom BLAST database location:

.. code-block:: bash

   # In bash
   export BLASTDB=/path/to/your/blastdb

   # Add to ~/.bashrc or ~/.bash_profile for persistence
   echo 'export BLASTDB=/path/to/your/blastdb' >> ~/.bashrc

Next Steps
----------

After installation, see the :doc:`quickstart` guide to begin using ePLACE.
