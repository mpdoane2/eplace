Phylogenetic Trees
==================

ePLACE uses IQTree to build phylogenetic trees from multiple sequence alignments.

Overview
--------

After aligning sequences, ePLACE builds phylogenetic trees to show the evolutionary relationships between query sequences and their taxonomic representatives. The tree building step:

* Uses maximum likelihood methods via IQTree
* Automatically selects the best-fit substitution model
* Performs ultrafast bootstrap analysis
* Labels tree tips with taxonomic information

IQTree Integration
------------------

ePLACE uses IQTree2 with automatic model selection and bootstrap support.

Model Selection
~~~~~~~~~~~~~~~

IQTree automatically selects the best-fit substitution model using ModelFinder:

.. code-block:: bash

   iqtree2 -s alignment.fasta -m MFP -B 1000 -T AUTO

Where:

* ``-m MFP``: ModelFinder Plus - tests and selects best model
* ``-B 1000``: 1000 ultrafast bootstrap replicates
* ``-T AUTO``: Automatic thread detection

Supported Models
~~~~~~~~~~~~~~~~

IQTree tests various nucleotide substitution models including:

* JC (Jukes-Cantor)
* F81
* K2P (Kimura 2-parameter)
* HKY
* TN (Tamura-Nei)
* TNe+I+G and variants
* GTR and variants

Using the API
-------------

Basic Tree Building
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.alignment import build_phylogenetic_tree

   # Build tree
   success = build_phylogenetic_tree(
       alignment_fasta=Path("aligned.fasta"),
       output_prefix=Path("tree"),
       num_threads=4
   )

   if success:
       print("Tree built successfully")
       print("Tree file: tree.treefile")
   else:
       print("Tree building failed")

Checking IQTree Availability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from eplace_lib.alignment import check_iqtree_available

   if check_iqtree_available():
       print("IQTree is available")
   else:
       print("IQTree not found - install IQTree to enable tree building")

Adding Taxonomic Labels
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.alignment import label_tree_with_taxonomy

   # Label tree tips with taxonomic information
   success = label_tree_with_taxonomy(
       tree_file=Path("tree.treefile"),
       output_tree=Path("tree_labeled.treefile"),
       blast_hits=filtered_hits,
       rank="genus"
   )

Tree Building in Workflows
---------------------------

Individual Workflow
~~~~~~~~~~~~~~~~~~~

In the individual workflow (``eplace search``):

1. Each query gets its own tree
2. Tree includes query + representative sequences
3. Tips are labeled with taxonomic information
4. Tree files are saved in query-specific directory

.. code-block:: bash

   eplace search query.fasta output_dir --tree-label-rank genus

Grouped Workflow
~~~~~~~~~~~~~~~~

In the grouped workflow (``eplace grouped``):

1. One tree per taxonomic group
2. Tree includes all queries in group + unique references
3. Shows relationships between multiple queries
4. Useful for comparative analysis

.. code-block:: bash

   eplace grouped query.fasta output_dir \
       --group-rank family \
       --tree-label-rank genus

Output Files
------------

IQTree produces several output files:

Primary Output
~~~~~~~~~~~~~~

* ``*.treefile`` - Best tree in Newick format (main output)
* ``*_labeled.treefile`` - Tree with taxonomic labels (ePLACE addition)

Supporting Files
~~~~~~~~~~~~~~~~

* ``*.iqtree`` - Full IQTree report with model selection and statistics
* ``*.log`` - Detailed log of tree building process
* ``*.bionj`` - Initial tree from BioNJ
* ``*.mldist`` - Maximum likelihood distance matrix
* ``*.model.gz`` - Model parameters (if applicable)
* ``*.splits.nex`` - Split support values in NEXUS format
* ``*.contree`` - Consensus tree (if bootstrap performed)
* ``*.ckp.gz`` - Checkpoint file (for resuming interrupted runs)

Tree File Format
~~~~~~~~~~~~~~~~

Trees are in Newick format:

.. code-block:: text

   (query_1:0.05,(ref_1:0.02,ref_2:0.03):0.04);

Labeled trees include taxonomic information:

.. code-block:: text

   (query_1:0.05,(ref_1|Escherichia:0.02,ref_2|Salmonella:0.03):0.04);

Visualizing Trees
-----------------

Using Python
~~~~~~~~~~~~

.. code-block:: python

   from Bio import Phylo
   import matplotlib.pyplot as plt

   # Read tree
   tree = Phylo.read("tree.treefile", "newick")

   # Draw tree
   fig = plt.figure(figsize=(10, 8))
   Phylo.draw(tree, do_show=False)
   plt.tight_layout()
   plt.savefig("tree.png", dpi=300)
   plt.show()

Using External Tools
~~~~~~~~~~~~~~~~~~~~

* **FigTree**: GUI application for viewing and annotating trees
* **iTOL**: Interactive Tree Of Life (web-based)
* **ggtree**: R package for tree visualization
* **ETE Toolkit**: Python framework for tree analysis and visualization

Example with ETE3:

.. code-block:: python

   from ete3 import Tree, TreeStyle

   # Read tree
   t = Tree("tree.treefile")

   # Style
   ts = TreeStyle()
   ts.show_leaf_name = True
   ts.show_branch_length = True
   ts.show_branch_support = True

   # Render
   t.render("tree.pdf", tree_style=ts)

Interpreting Trees
------------------

Branch Lengths
~~~~~~~~~~~~~~

* Represent evolutionary distance (substitutions per site)
* Longer branches = more evolutionary change
* Scale bar shows units

Bootstrap Support
~~~~~~~~~~~~~~~~~

* Numbers at nodes indicate support (0-100)
* >95: Strong support
* 70-95: Moderate support
* <70: Weak support

Tree Topology
~~~~~~~~~~~~~

* Sister taxa are more closely related
* Deeper nodes = older divergence
* Monophyletic groups share common ancestor

Troubleshooting
---------------

IQTree not found
~~~~~~~~~~~~~~~~

If you get "IQTree is not available":

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install iqtree

   # macOS
   brew install iqtree

   # Conda
   conda install -c bioconda iqtree

Tree building fails
~~~~~~~~~~~~~~~~~~~

Common causes:

1. **Insufficient sequences**: Need ≥3 sequences for tree
2. **Poor alignment**: Check alignment quality first
3. **Identical sequences**: Remove duplicates
4. **No variation**: All sequences too similar

Tree building too slow
~~~~~~~~~~~~~~~~~~~~~~

For faster tree building:

1. Increase ``--num-threads``
2. Reduce bootstrap replicates (not recommended for publication)
3. Use simpler models
4. Reduce number of sequences

Strange tree topology
~~~~~~~~~~~~~~~~~~~~~

If tree structure seems incorrect:

1. Check alignment quality
2. Verify sequences are homologous
3. Check for contamination or misidentification
4. Consider longer sequences for better resolution
5. Try different taxonomic ranks for representatives

Advanced Usage
--------------

Custom IQTree Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

Modify tree building parameters:

.. code-block:: python

   import subprocess
   from pathlib import Path

   def custom_tree_build(alignment: Path, prefix: Path):
       """Build tree with custom IQTree parameters."""
       cmd = [
           "iqtree2",
           "-s", str(alignment),
           "-pre", str(prefix),
           "-m", "GTR+I+G",  # Specific model
           "-B", "2000",      # More bootstrap replicates
           "-alrt", "1000",   # SH-aLRT test
           "-T", "8"
       ]
       subprocess.run(cmd, check=True)

Parsing Tree Files
~~~~~~~~~~~~~~~~~~

Extract information from trees:

.. code-block:: python

   from Bio import Phylo

   # Read tree
   tree = Phylo.read("tree.treefile", "newick")

   # Get all tips
   tips = tree.get_terminals()
   print(f"Number of tips: {len(tips)}")

   # Calculate tree height
   height = tree.total_branch_length()
   print(f"Tree height: {height:.4f}")

   # Find specific clade
   for clade in tree.find_clades():
       if clade.name and "query" in clade.name:
           print(f"Found query: {clade.name}")

Tree Comparison
~~~~~~~~~~~~~~~

Compare multiple trees:

.. code-block:: python

   from Bio import Phylo
   from Bio.Phylo.Consensus import majority_consensus

   # Read multiple trees
   trees = list(Phylo.parse("trees.nexus", "nexus"))

   # Build consensus tree
   consensus = majority_consensus(trees, cutoff=0.5)

   # Write consensus
   Phylo.write(consensus, "consensus.tree", "newick")

Rerooting Trees
~~~~~~~~~~~~~~~

Change tree root:

.. code-block:: python

   from Bio import Phylo

   # Read tree
   tree = Phylo.read("tree.treefile", "newick")

   # Reroot on outgroup
   outgroup = tree.find_any(name="outgroup_name")
   tree.root_with_outgroup(outgroup)

   # Write rerooted tree
   Phylo.write(tree, "rerooted_tree.treefile", "newick")

Best Practices
--------------

1. **Alignment Quality**: Always check alignment before tree building
2. **Bootstrap Support**: Use bootstrap to assess confidence
3. **Model Selection**: Let IQTree select best model automatically
4. **Outgroups**: Include outgroup if possible for rooting
5. **Visualization**: View trees to catch obvious errors
6. **Documentation**: Record all parameters for reproducibility

Performance Considerations
--------------------------

Tree building time depends on:

* Number of sequences
* Sequence length
* Model complexity
* Number of bootstrap replicates
* Number of threads
* CPU speed

Typical timings:

* 10 sequences × 500bp: ~10-30 seconds
* 50 sequences × 1000bp: ~1-5 minutes
* 100 sequences × 2000bp: ~10-30 minutes

Using maximum threads significantly improves speed.

Statistical Considerations
--------------------------

Branch Support
~~~~~~~~~~~~~~

* **UFBoot**: Ultrafast bootstrap approximation (default)
* **Standard Bootstrap**: Classic but slower
* **SH-aLRT**: Shimodaira-Hasegawa-like approximate likelihood ratio test

All methods assess confidence in tree topology.

Model Selection Criteria
~~~~~~~~~~~~~~~~~~~~~~~~

IQTree uses:

* **BIC**: Bayesian Information Criterion (default)
* **AIC**: Akaike Information Criterion
* **AICc**: Corrected AIC

BIC generally preferred to avoid over-parameterization.

See Also
--------

* :doc:`alignment` - Sequence alignment documentation
* :doc:`workflows` - Complete workflow documentation
* :doc:`blast_workflow` - Full pipeline from BLAST to trees
* `IQTree documentation <http://www.iqtree.org/doc/>`_
* `Phylogenetics guide <http://www.phylo.org/>`_
