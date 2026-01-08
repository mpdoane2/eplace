API Reference
=============

This page documents the ePLACE Python API for programmatic access.

Core Modules
------------

eplace_lib.blast_analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: eplace_lib.blast_analysis
   :members:
   :undoc-members:
   :show-inheritance:

eplace_lib.taxonomy
~~~~~~~~~~~~~~~~~~~

.. automodule:: eplace_lib.taxonomy
   :members:
   :undoc-members:
   :show-inheritance:

eplace_lib.sequences
~~~~~~~~~~~~~~~~~~~~

.. automodule:: eplace_lib.sequences
   :members:
   :undoc-members:
   :show-inheritance:

eplace_lib.alignment
~~~~~~~~~~~~~~~~~~~~

.. automodule:: eplace_lib.alignment
   :members:
   :undoc-members:
   :show-inheritance:

eplace_lib.ncbi_download
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: eplace_lib.ncbi_download
   :members:
   :undoc-members:
   :show-inheritance:

eplace_lib.cli
~~~~~~~~~~~~~~

.. automodule:: eplace_lib.cli
   :members:
   :undoc-members:
   :show-inheritance:

Quick Examples
--------------

BLAST Analysis
~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib import run_blast_search, process_blast_results_for_taxonomy

   # Run BLAST search with filtering
   success, filtered_hits = run_blast_search(
       query_fasta=Path("query.fasta"),
       output_file=Path("blast_results.txt"),
       min_identity=90.0,
       min_coverage=80.0
   )

   # Extract representative sequences
   results = process_blast_results_for_taxonomy(
       blast_hits=filtered_hits,
       output_dir=Path("output"),
       rank="genus"
   )

Database Download
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from eplace_lib import setup_ncbi_database

   # Download the core_nt database
   success, message = setup_ncbi_database()
   print(f"Success: {success}, Message: {message}")

FASTA Reading
~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.blast_analysis import FastaReader

   # Read sequences
   sequences = FastaReader.read_fasta(Path("input.fasta"))

   # Get sequence lengths
   lengths = FastaReader.get_sequence_lengths(Path("input.fasta"))

Sequence Alignment
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.alignment import align_sequences, build_phylogenetic_tree

   # Align sequences
   success = align_sequences(
       input_fasta=Path("sequences.fasta"),
       output_fasta=Path("aligned.fasta"),
       num_threads=4
   )

   # Build tree
   success = build_phylogenetic_tree(
       alignment_fasta=Path("aligned.fasta"),
       output_prefix=Path("tree"),
       num_threads=4
   )

Data Structures
---------------

BlastHit
~~~~~~~~

Represents a single BLAST hit with the following attributes:

* ``query_id``: Query sequence identifier
* ``subject_id``: Subject (database) sequence identifier
* ``percent_identity``: Percentage of identical matches
* ``alignment_length``: Length of alignment
* ``query_length``: Length of query sequence
* ``subject_length``: Length of subject sequence
* ``query_start``: Start position in query
* ``query_end``: End position in query
* ``subject_start``: Start position in subject
* ``subject_end``: End position in subject
* ``evalue``: Expectation value
* ``bit_score``: Bit score
* ``query_coverage``: Percentage of query covered by alignment
* ``subject_taxonomy``: Dictionary containing taxonomic information (phylum, class, order, family, genus, species)

Example usage:

.. code-block:: python

   from eplace_lib.blast_analysis import BlastHit

   # Create a BlastHit
   hit = BlastHit(
       query_id="query1",
       subject_id="NC_001234.5",
       percent_identity=95.5,
       alignment_length=500,
       query_length=550,
       subject_length=5000,
       query_start=1,
       query_end=500,
       subject_start=100,
       subject_end=599,
       evalue=1e-100,
       bit_score=900,
       query_coverage=90.9,
       subject_taxonomy={"genus": "Escherichia", "species": "coli"}
   )

Common Workflows
----------------

Complete BLAST to Tree Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib import (
       run_blast_search,
       process_blast_results_for_taxonomy,
   )
   from eplace_lib.sequences import trim_sequences_to_blast_coordinates
   from eplace_lib.alignment import align_sequences, build_phylogenetic_tree

   # Step 1: BLAST search
   success, filtered_hits = run_blast_search(
       query_fasta=Path("query.fasta"),
       output_file=Path("blast_results.txt"),
       min_identity=90.0,
       min_coverage=80.0,
       num_threads=4
   )

   # Step 2: Extract representatives
   results = process_blast_results_for_taxonomy(
       blast_hits=filtered_hits,
       output_dir=Path("output"),
       rank="genus"
   )

   # Step 3: Process each query
   for query_id, fasta_path in results.items():
       # Trim sequences
       trimmed_path = fasta_path.parent / f"{query_id}_trimmed.fasta"
       trim_sequences_to_blast_coordinates(
           input_fasta=fasta_path,
           output_fasta=trimmed_path,
           blast_hits=filtered_hits
       )
       
       # Align sequences
       aligned_path = fasta_path.parent / f"{query_id}_aligned.fasta"
       align_sequences(
           input_fasta=trimmed_path,
           output_fasta=aligned_path,
           num_threads=4
       )
       
       # Build tree
       tree_prefix = fasta_path.parent / f"{query_id}_tree"
       build_phylogenetic_tree(
           alignment_fasta=aligned_path,
           output_prefix=tree_prefix,
           num_threads=4
       )

Custom BLAST Parameters
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.blast_analysis import BlastRunner

   runner = BlastRunner()

   # Run BLAST with custom parameters
   success = runner.run_blastn(
       query_fasta=Path("query.fasta"),
       output_file=Path("blast_results.txt"),
       database="core_nt",
       num_threads=8,
       max_target_seqs=500,
       evalue=1e-10,
       word_size=11
   )

   # Parse and filter results
   hits = runner.parse_blast_results(Path("blast_results.txt"))
   filtered_hits = runner.filter_blast_hits(
       hits,
       min_identity=95.0,
       min_coverage=90.0
   )

Working with Taxonomic Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from eplace_lib.taxonomy import TaxonomyExtractor

   extractor = TaxonomyExtractor()

   # Group hits by query
   grouped_hits = extractor.group_hits_by_query(blast_hits)

   # Select representatives at different ranks
   for query_id, query_hits in grouped_hits.items():
       # At genus level
       genus_reps = extractor.select_representatives_by_rank(
           hits=query_hits,
           rank="genus",
           max_per_rank=1
       )
       
       # At species level
       species_reps = extractor.select_representatives_by_rank(
           hits=query_hits,
           rank="species",
           max_per_rank=2
       )

Error Handling
--------------

Most functions return success indicators and provide error messages:

.. code-block:: python

   from pathlib import Path
   from eplace_lib import run_blast_search

   success, result = run_blast_search(
       query_fasta=Path("query.fasta"),
       output_file=Path("output.txt"),
       min_identity=90.0,
       min_coverage=80.0
   )

   if not success:
       print(f"BLAST failed: {result}")
   else:
       print(f"Found {len(result)} hits")

For functions that don't return tuples, check return values:

.. code-block:: python

   from pathlib import Path
   from eplace_lib.alignment import align_sequences

   success = align_sequences(
       input_fasta=Path("sequences.fasta"),
       output_fasta=Path("aligned.fasta")
   )

   if not success:
       print("Alignment failed")

Type Hints
----------

ePLACE uses type hints throughout the codebase for better IDE support:

.. code-block:: python

   from pathlib import Path
   from typing import List, Dict, Tuple
   from eplace_lib.blast_analysis import BlastHit

   def process_hits(
       hits: List[BlastHit],
       min_identity: float = 90.0
   ) -> Tuple[bool, List[BlastHit]]:
       """Process BLAST hits with type hints."""
       filtered = [h for h in hits if h.percent_identity >= min_identity]
       return True, filtered

Logging
-------

ePLACE uses Python's logging module. Configure logging in your scripts:

.. code-block:: python

   import logging

   # Configure logging
   logging.basicConfig(
       level=logging.INFO,
       format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
   )

   # Now run ePLACE functions
   from eplace_lib import run_blast_search

Advanced Usage
--------------

Custom Database Management
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from eplace_lib.ncbi_download import NCBIDownloader

   downloader = NCBIDownloader()

   # Get database directory
   db_dir = downloader.get_blastdb_directory()

   # Check if database exists
   exists = downloader.check_database_exists()

   # Get available files
   files = downloader.get_available_files()

   # Download specific file
   downloader.download_file('core_nt.00.tar.gz', db_dir)

Sequence Extraction
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.taxonomy import SequenceExtractor

   extractor = SequenceExtractor()

   # Extract specific sequences
   success = extractor.extract_sequences(
       sequence_ids=["NC_001234.5", "NC_005678.9"],
       output_fasta=Path("extracted.fasta"),
       database="core_nt"
   )

See Also
--------

* :doc:`quickstart` - Quick start guide with examples
* :doc:`workflows` - Workflow documentation
* :doc:`blast_workflow` - Detailed BLAST workflow guide
