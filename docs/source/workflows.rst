Workflows
=========

ePLACE provides two main workflow types for analyzing environmental DNA sequences: **Individual** and **Grouped** workflows.

Workflow Overview
-----------------

Both workflows follow a similar pipeline but differ in how they organize and analyze query sequences.

Common Pipeline Steps
~~~~~~~~~~~~~~~~~~~~~

1. **BLAST Search** - Search query sequences against NCBI database
2. **Filter Results** - Apply identity and coverage thresholds
3. **Taxonomic Classification** - Classify hits using TaxonKit
4. **Representative Selection** - Select representative sequences per taxonomic rank
5. **Sequence Extraction** - Extract sequences from BLAST database
6. **Sequence Trimming** - Trim to aligned regions based on BLAST coordinates
7. **Multiple Sequence Alignment** - Align using MAFFT (optional)
8. **Phylogenetic Tree Building** - Build trees using IQTree (optional)

Individual Workflow
-------------------

The individual workflow (``eplace blast``) processes each query sequence independently.

When to Use
~~~~~~~~~~~

* Analyzing diverse sequences that may not be related
* Need independent phylogenetic context for each sequence
* Want clear, per-sequence results
* Queries span multiple taxonomic groups

Process
~~~~~~~

.. code-block:: text

   For each query sequence:
   1. BLAST search → filter hits
   2. Classify hits taxonomically
   3. Select representatives at specified rank
   4. Create query-specific directory
   5. Extract and trim sequences
   6. Align sequences (query + representatives)
   7. Build phylogenetic tree

Output Structure
~~~~~~~~~~~~~~~~

.. code-block:: text

   output_dir/
   ├── blast_results.txt
   ├── blast_results_annotated.txt
   ├── query1/
   │   ├── query1_representatives.fasta
   │   ├── query1_with_query.fasta
   │   ├── query1_trimmed.fasta
   │   ├── query1_aligned.fasta
   │   └── query1_tree.treefile
   ├── query2/
   │   └── ...
   └── ...

Usage Example
~~~~~~~~~~~~~

.. code-block:: bash

   # Basic individual workflow
   eplace blast queries.fasta output_dir

   # With custom parameters
   eplace blast queries.fasta output_dir \
       --rank genus \
       --min-identity 95 \
       --min-coverage 85 \
       --num-threads 4

   # BLAST only, no alignment
   eplace blast queries.fasta output_dir --skip-alignment

Grouped Workflow
----------------

The grouped workflow (``eplace grouped``) combines queries by taxonomic classification for joint analysis.

When to Use
~~~~~~~~~~~

* Analyzing related sequences from similar taxonomic groups
* Want to see multiple queries in single phylogenetic context
* Need comparative analysis of related sequences
* Want to reduce computational overhead

Process
~~~~~~~

.. code-block:: text

   1. BLAST search all queries → filter hits
   2. Classify hits taxonomically
   3. Select representatives for each query at specified rank
   4. Group queries by taxonomic rank (e.g., class, order, family)
   5. For each group:
      a. Combine all queries in group
      b. Collect unique reference sequences
      c. Remove redundant references
      d. Verify alignment consistency
      e. Extract and trim sequences
      f. Align all sequences together
      g. Build single phylogenetic tree

Output Structure
~~~~~~~~~~~~~~~~

.. code-block:: text

   output_dir/
   ├── blast_results.txt
   ├── blast_results_annotated.txt
   ├── query1/
   │   └── query1_representatives.fasta
   ├── query2/
   │   └── query2_representatives.fasta
   ├── Bacteria_Proteobacteria/
   │   ├── Bacteria_Proteobacteria_combined.fasta
   │   ├── Bacteria_Proteobacteria_trimmed.fasta
   │   ├── Bacteria_Proteobacteria_aligned.fasta
   │   └── Bacteria_Proteobacteria_tree.treefile
   ├── Bacteria_Firmicutes/
   │   └── ...
   └── ...

Usage Example
~~~~~~~~~~~~~

.. code-block:: bash

   # Basic grouped workflow (groups by class)
   eplace grouped queries.fasta output_dir

   # Group by different rank
   eplace grouped queries.fasta output_dir --group-rank order

   # Specify both representative and grouping ranks
   eplace grouped queries.fasta output_dir \
       --rank genus \
       --group-rank family

   # With custom parameters
   eplace grouped queries.fasta output_dir \
       --rank genus \
       --group-rank family \
       --min-identity 95 \
       --min-coverage 85 \
       --num-threads 4

Comparison
----------

+-----------------------+-------------------------+-------------------------+
| Feature               | Individual              | Grouped                 |
+=======================+=========================+=========================+
| Analysis unit         | Per query               | Per taxonomic group     |
+-----------------------+-------------------------+-------------------------+
| Trees generated       | One per query           | One per group           |
+-----------------------+-------------------------+-------------------------+
| Alignment scope       | Query + its refs        | All queries + unique    |
|                       |                         | refs in group           |
+-----------------------+-------------------------+-------------------------+
| Best for              | Diverse sequences       | Related sequences       |
+-----------------------+-------------------------+-------------------------+
| Computational cost    | Higher (more trees)     | Lower (fewer trees)     |
+-----------------------+-------------------------+-------------------------+
| Interpretation        | Independent per query   | Comparative across      |
|                       |                         | queries                 |
+-----------------------+-------------------------+-------------------------+

Workflow Parameters
-------------------

Taxonomic Rank Selection
~~~~~~~~~~~~~~~~~~~~~~~~~

Both workflows use taxonomic ranks to organize results:

* ``--rank``: Level at which to select representative sequences

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: ``genus``
  * Higher ranks = fewer, more diverse representatives
  * Lower ranks = more, more specific representatives

* ``--group-rank`` (grouped only): Level at which to group queries

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: ``class``
  * Determines how queries are combined
  * Should typically be equal to or higher than ``--rank``

* ``--tree-label-rank``: Level at which to label tree tips

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: ``genus``
  * Determines taxonomic labels on phylogenetic tree

Filtering Parameters
~~~~~~~~~~~~~~~~~~~~

Control which BLAST hits are included:

* ``--min-identity``: Minimum percent identity (default: 90.0)

  * Range: 0-100
  * Higher = more stringent matching
  * Typical values: 90-99

* ``--min-coverage``: Minimum query coverage (default: 80.0)

  * Range: 0-100
  * Percentage of query sequence covered by alignment
  * Higher = require longer alignments

Database Parameters
~~~~~~~~~~~~~~~~~~~

Configure BLAST database usage:

* ``--database``: BLAST database name (default: core_nt)

  * Standard databases: nt, core_nt, refseq_rna, etc.
  * Must be installed in BLASTDB path

* ``--blastdb-path``: Custom path to BLAST databases

  * Overrides $BLASTDB environment variable
  * Use for non-standard locations

Performance Parameters
~~~~~~~~~~~~~~~~~~~~~~

Optimize analysis speed:

* ``--num-threads``: Number of CPU threads to use (default: 1)

  * Used by BLAST, MAFFT, and IQTree
  * Set to number of available cores for speed
  * Higher values = faster but more resource intensive

Workflow Control
~~~~~~~~~~~~~~~~

Fine-tune workflow behavior:

* ``--skip-alignment``: Skip MAFFT alignment and IQTree tree building

  * Useful for BLAST-only analysis
  * Saves time and computational resources
  * Still performs extraction and trimming

* ``--overwrite-existing-blast``: Re-run BLAST even if results exist

  * By default, existing BLAST results are reused
  * Use when BLAST parameters change

* ``--alignment-tolerance`` (grouped only): Maximum coordinate difference for alignment consistency (default: 50)

  * Checks if BLAST hits align to similar regions
  * Larger values = more permissive grouping
  * Smaller values = stricter alignment requirements

* ``--output-classification``: Path to output classification TSV file

  * Generates summary table of taxonomic classifications
  * Useful for downstream analysis

Advanced Usage
--------------

Chaining Analyses
~~~~~~~~~~~~~~~~~

You can run different analyses on the same BLAST results:

.. code-block:: bash

   # First run: BLAST only
   eplace blast queries.fasta output1 --skip-alignment

   # Analyze at different ranks using same BLAST results
   eplace blast queries.fasta output2 --rank genus
   eplace blast queries.fasta output3 --rank species

Batch Processing
~~~~~~~~~~~~~~~~

Process multiple query files:

.. code-block:: bash

   for file in queries/*.fasta; do
       base=$(basename "$file" .fasta)
       eplace blast "$file" "output_$base" --num-threads 8
   done

Custom Filtering Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~

Apply different stringency levels:

.. code-block:: bash

   # High stringency
   eplace blast queries.fasta output_strict \
       --min-identity 98 \
       --min-coverage 95

   # Medium stringency
   eplace blast queries.fasta output_medium \
       --min-identity 95 \
       --min-coverage 85

   # Low stringency
   eplace blast queries.fasta output_relaxed \
       --min-identity 90 \
       --min-coverage 75

Hierarchical Grouping
~~~~~~~~~~~~~~~~~~~~~~

Analyze at multiple taxonomic levels:

.. code-block:: bash

   # Group by class
   eplace grouped queries.fasta output_class --group-rank class

   # Group by order (more specific)
   eplace grouped queries.fasta output_order --group-rank order

   # Group by family (even more specific)
   eplace grouped queries.fasta output_family --group-rank family

Best Practices
--------------

Choosing Parameters
~~~~~~~~~~~~~~~~~~~

1. **Start with defaults** to get a baseline
2. **Adjust identity/coverage** based on sequence quality
3. **Choose rank based on analysis goals**:

   * Phylum/Class: Broad overview
   * Order/Family: Balanced detail
   * Genus/Species: Specific identification

4. **Use grouped workflow** when queries are known to be related
5. **Use ``--num-threads``** to speed up analysis

Quality Control
~~~~~~~~~~~~~~~

1. Check BLAST results for reasonable number of hits
2. Verify taxonomic classifications make sense
3. Inspect alignments for quality
4. Review phylogenetic trees for expected topology

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

1. Use ``--skip-alignment`` for exploratory analysis
2. Increase ``--num-threads`` on multi-core systems
3. Use grouped workflow to reduce number of trees
4. Filter query sequences before analysis
5. Consider splitting very large query files

Troubleshooting
---------------

No sequences in group
~~~~~~~~~~~~~~~~~~~~~

If grouped workflow creates empty groups:

* Queries may be too diverse
* Try lower ``--group-rank`` (e.g., phylum instead of class)
* Check taxonomic classifications of queries

Alignment consistency errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If grouped workflow reports alignment inconsistencies:

* BLAST hits align to different regions of references
* Increase ``--alignment-tolerance``
* Or use individual workflow instead

Tree building fails
~~~~~~~~~~~~~~~~~~~

If IQTree fails to build tree:

* Check alignment quality
* Ensure sufficient sequences (≥3)
* Verify sequences have overlapping regions
* Try ``--skip-alignment`` to see raw sequences

See Also
--------

* :doc:`cli` - Complete command-line reference
* :doc:`blast_workflow` - Detailed BLAST workflow documentation
* :doc:`api` - Python API for custom workflows
