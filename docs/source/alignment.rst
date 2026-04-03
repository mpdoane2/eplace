Sequence Alignment
==================

ePLACE uses MAFFT for multiple sequence alignment of query sequences and their taxonomic representatives.

Overview
--------

After extracting representative sequences from BLAST results, ePLACE aligns them together with the query sequence to prepare for phylogenetic analysis. The alignment step:

* Combines query and reference sequences
* Uses MAFFT with auto-orientation (``--adjustdirection``)
* Handles sequences in different orientations
* Produces aligned FASTA files suitable for tree building

MAFFT Integration
-----------------

ePLACE uses MAFFT's automatic alignment mode with several key features:

Auto-orientation
~~~~~~~~~~~~~~~~

The ``--adjustdirection`` flag automatically detects and corrects reverse-complemented sequences:

.. code-block:: bash

   mafft --auto --adjustdirection --thread N input.fasta > output.fasta

This is crucial for environmental DNA sequences which may be in either orientation.

Automatic Algorithm Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MAFFT's ``--auto`` flag selects the most appropriate algorithm based on:

* Number of sequences
* Sequence lengths
* Available memory

This provides a good balance of speed and accuracy for most datasets.

Using the API
-------------

Basic Alignment
~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.alignment import align_sequences

   # Align sequences
   success = align_sequences(
       input_fasta=Path("sequences.fasta"),
       output_fasta=Path("aligned.fasta"),
       num_threads=4
   )

   if success:
       print("Alignment successful")
   else:
       print("Alignment failed")

Checking MAFFT Availability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from eplace_lib.alignment import check_mafft_available

   if check_mafft_available():
       print("MAFFT is available")
   else:
       print("MAFFT not found - install MAFFT to enable alignment")

Alignment in Workflows
----------------------

Individual Workflow
~~~~~~~~~~~~~~~~~~~

In the individual workflow (``eplace search``):

1. Query sequence is combined with its representatives
2. Combined FASTA is trimmed to aligned regions
3. Trimmed sequences are aligned with MAFFT
4. Aligned sequences are used for tree building

.. code-block:: bash

   eplace search query.fasta output_dir --num-threads 4

Grouped Workflow
~~~~~~~~~~~~~~~~

In the grouped workflow (``eplace grouped``):

1. Multiple queries and their unique references are combined
2. Combined FASTA is trimmed to aligned regions
3. All sequences in the group are aligned together
4. Aligned sequences are used for phylogenetic tree

.. code-block:: bash

   eplace grouped query.fasta output_dir --num-threads 4

Skipping Alignment
------------------

If you only need BLAST results without alignment:

.. code-block:: bash

   eplace search query.fasta output_dir --skip-alignment

This will:

* Perform BLAST search
* Extract representative sequences
* Skip MAFFT alignment
* Skip tree building

Output Files
------------

Alignment produces the following files:

* ``*_aligned.fasta`` - Multiple sequence alignment in FASTA format
* Sequences are in the same order as input
* Gap characters (``-``) indicate alignment positions
* Reversed sequences are marked with ``_R_`` prefix (MAFFT convention)

Example aligned FASTA:

.. code-block:: text

   >query_sequence_1
   ATGC-ATGCATGC
   >representative_1
   ATGCGATGCATGC
   >_R_representative_2
   ATGC-ATGCATGC

Tree Labeling
~~~~~~~~~~~~~

When phylogenetic trees are labeled with taxonomic names, sequences that were reverse-complemented by MAFFT (marked with ``_R_`` prefix) will have ``_R`` appended to their taxonomic label. This allows you to identify which sequences were reverse-complemented during alignment.

Example tree before relabeling:

.. code-block:: text

   (MZ387488.1:0.1,_R_CP123456.1:0.2,query:0.0);

Example tree after relabeling:

.. code-block:: text

   (Salmonella:0.1,Escherichia_R:0.2,query:0.0);

The ``_R`` suffix indicates that the Escherichia sequence was reverse-complemented during alignment.

Troubleshooting
---------------

MAFFT not found
~~~~~~~~~~~~~~~

If you get "MAFFT is not available":

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install mafft

   # macOS
   brew install mafft

   # Conda
   conda install -c bioconda mafft

Alignment takes too long
~~~~~~~~~~~~~~~~~~~~~~~~~

For large datasets:

1. Increase ``--num-threads`` for parallelization
2. Consider filtering to fewer representative sequences
3. Use grouped workflow to reduce number of alignments
4. Pre-filter sequences by length or quality

Out of memory
~~~~~~~~~~~~~

If MAFFT runs out of memory:

1. Process fewer sequences at a time
2. Reduce number of representatives per rank
3. Use a machine with more RAM
4. Consider sequence length limits

Poor alignment quality
~~~~~~~~~~~~~~~~~~~~~~

If alignments appear poor:

1. Check sequence quality and length
2. Ensure sequences are from related organisms
3. Verify BLAST filtering parameters are appropriate
4. Consider manual curation of representatives

Advanced Usage
--------------

Custom MAFFT Parameters
~~~~~~~~~~~~~~~~~~~~~~~

While ePLACE uses sensible defaults, you can modify the alignment code to use custom MAFFT parameters:

.. code-block:: python

   import subprocess
   from pathlib import Path

   def custom_align(input_fasta: Path, output_fasta: Path):
       """Custom alignment with specific MAFFT parameters."""
       cmd = [
           "mafft",
           "--maxiterate", "1000",
           "--localpair",  # L-INS-i algorithm
           "--adjustdirection",
           "--thread", "8",
           str(input_fasta)
       ]
       
       with open(output_fasta, 'w') as f:
           subprocess.run(cmd, stdout=f, check=True)

Alignment Quality Assessment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assess alignment quality:

.. code-block:: python

   from Bio import AlignIO

   # Read alignment
   alignment = AlignIO.read("aligned.fasta", "fasta")

   # Calculate statistics
   num_seqs = len(alignment)
   align_length = alignment.get_alignment_length()
   
   # Count gaps
   gap_counts = [seq.seq.count('-') for seq in alignment]
   avg_gaps = sum(gap_counts) / num_seqs
   
   print(f"Sequences: {num_seqs}")
   print(f"Alignment length: {align_length}")
   print(f"Average gaps: {avg_gaps:.1f}")

Post-processing Alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trim poorly aligned regions:

.. code-block:: python

   from Bio import AlignIO
   from Bio.Align import MultipleSeqAlignment

   # Read alignment
   alignment = AlignIO.read("aligned.fasta", "fasta")

   # Trim columns with >50% gaps
   trimmed_alignment = []
   for i in range(alignment.get_alignment_length()):
       column = alignment[:, i]
       gap_fraction = column.count('-') / len(column)
       if gap_fraction < 0.5:
           trimmed_alignment.append(column)

   # Save trimmed alignment
   # (implementation depends on your needs)

Best Practices
--------------

1. **Quality Control**: Check input sequences before alignment
2. **Threading**: Use multiple threads (``--num-threads``) for speed
3. **Memory**: Monitor memory usage for large alignments
4. **Validation**: Visually inspect alignments when possible
5. **Documentation**: Record alignment parameters for reproducibility

Performance Considerations
--------------------------

Alignment speed depends on:

* Number of sequences
* Sequence lengths
* Algorithm selected by MAFFT
* Number of threads
* Available memory

Typical timings:

* 10 sequences × 500bp: < 1 second
* 50 sequences × 1000bp: ~5-10 seconds
* 100 sequences × 2000bp: ~30-60 seconds

See Also
--------

* :doc:`phylogenetic_trees` - Tree building from alignments
* :doc:`workflows` - Complete workflow documentation
* :doc:`blast_workflow` - BLAST to alignment pipeline
* `MAFFT documentation <https://mafft.cbrc.jp/alignment/software/>`_
