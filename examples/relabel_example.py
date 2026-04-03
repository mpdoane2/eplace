#!/usr/bin/env python
"""
Example script demonstrating the use of the 'eplace relabel' command.

This example shows how to relabel a phylogenetic tree with taxonomic names
from BLAST results.

Requirements:
- BLAST results file (tabular format with taxonomy information)
- Phylogenetic tree file (Newick format)
- pytaxonkit installed

Usage:
    python relabel_example.py
    
    # Or use the CLI directly:
    eplace relabel blast_results.txt input.treefile output.treefile --rank genus
"""

import subprocess
from pathlib import Path


def create_sample_data():
    """
    Create sample BLAST and tree files for demonstration.
    
    Note: In real usage, these files would come from actual BLAST searches
    and phylogenetic tree building tools.
    """
    
    # Create sample BLAST output (simplified format)
    # Format: query_id, subject_id, pident, length, qlen, slen, qstart, qend,
    #         sstart, send, evalue, bitscore, staxid, staxids
    blast_data = """query1\tNC_000913.3\t98.5\t500\t1000\t4641652\t1\t500\t1000\t1500\t0.0\t925\t562\t562
query1\tNC_002695.2\t95.2\t480\t1000\t5498578\t1\t480\t2000\t2480\t1e-180\t850\t562\t562
query1\tNC_004337.2\t93.8\t475\t1000\t4857432\t1\t475\t3000\t3475\t1e-175\t820\t216592\t216592
query2\tNC_003197.2\t97.5\t520\t1100\t5231428\t1\t520\t4000\t4520\t0.0\t950\t224308\t224308
query2\tNC_000964.3\t94.1\t510\t1100\t4215606\t1\t510\t5000\t5510\t1e-185\t880\t224308\t224308
"""
    
    # Create sample tree (Newick format)
    # Tree contains the accessions from the BLAST results
    tree_data = """((NC_000913.3:0.05,NC_002695.2:0.08):0.1,(NC_004337.2:0.12,(NC_003197.2:0.06,NC_000964.3:0.09):0.15):0.2);"""
    
    # Write files
    blast_file = Path("sample_blast_results.txt")
    tree_file = Path("sample_tree.treefile")
    
    with open(blast_file, 'w') as f:
        f.write(blast_data)
    
    with open(tree_file, 'w') as f:
        f.write(tree_data)
    
    print("✓ Created sample data files:")
    print(f"  - {blast_file}")
    print(f"  - {tree_file}")
    print()
    
    return blast_file, tree_file


def run_relabel_example(blast_file, tree_file, rank='genus'):
    """
    Run the eplace relabel command with sample data.
    
    Args:
        blast_file: Path to BLAST results file
        tree_file: Path to input tree file
        rank: Taxonomic rank for labeling (default: genus)
    """
    output_tree = Path(f"relabeled_tree_{rank}.treefile")
    
    print(f"Running relabel command with rank: {rank}")
    print("=" * 60)
    
    # Build the command
    cmd = [
        'eplace', 'relabel',
        str(blast_file),
        str(tree_file),
        str(output_tree),
        '--rank', rank
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print()
    
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            print("✓ Relabeling completed successfully!")
            print(f"✓ Output tree: {output_tree}")
            print()
            
            # Show the relabeled tree
            if output_tree.exists():
                with open(output_tree, 'r') as f:
                    tree_content = f.read()
                print("Relabeled tree content:")
                print(tree_content)
                print()
        else:
            print(f"✗ Error: {result.stderr}")
            
    except subprocess.TimeoutExpired:
        print("✗ Command timed out")
    except FileNotFoundError:
        print("✗ eplace command not found. Make sure it's installed:")
        print("   pip install -e .")


def demonstrate_all_ranks():
    """
    Demonstrate relabeling with different taxonomic ranks.
    """
    print("=" * 70)
    print("ePLACE Relabel Command Example")
    print("=" * 70)
    print()
    
    # Create sample data
    blast_file, tree_file = create_sample_data()
    
    # Demonstrate different ranks
    ranks = ['genus', 'species', 'family']
    
    for rank in ranks:
        print(f"\n{'=' * 70}")
        print(f"Example: Relabeling with {rank} rank")
        print(f"{'=' * 70}\n")
        
        if rank == 'species':
            print("Note: Species rank uses 'genus species' format for clarity")
            print()
        
        run_relabel_example(blast_file, tree_file, rank)
    
    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)


def show_usage_tips():
    """
    Display usage tips and best practices.
    """
    print("\n" + "=" * 70)
    print("Usage Tips:")
    print("=" * 70)
    print("""
1. BLAST Output Format:
   - Must be in tabular format with taxonomy information
   - Should include columns: qseqid, sseqid, pident, length, qlen, slen,
     qstart, qend, sstart, send, evalue, bitscore, staxid, staxids
   
2. Tree File Format:
   - Must be in Newick format
   - Sequence IDs in the tree should match the subject IDs from BLAST results
   
3. Taxonomic Ranks:
   - Available ranks: phylum, class, order, family, genus, species
   - Species rank automatically uses "genus species" format for clarity
   
4. Common Use Cases:
   - Genus labeling: Good balance between specificity and readability
   - Species labeling: Most specific, uses binomial nomenclature
   - Family/Order labeling: Useful for higher-level groupings
   
5. Command Examples:
   # Basic usage with genus
   eplace relabel blast.txt tree.nwk output.nwk --rank genus
   
   # Species labeling (genus + species)
   eplace relabel blast.txt tree.nwk output.nwk --rank species
   
   # Family labeling
   eplace relabel blast.txt tree.nwk output.nwk --rank family

6. Integration with Other Commands:
   # First run search workflow
   eplace search query.fasta output_dir --rank genus
   
   # Then relabel the tree with different rank
   # Use blast_results.txt (BLAST) or mmseqs_results.txt (MMseqs2)
   eplace relabel output_dir/blast_results.txt \\
          output_dir/query/tree.treefile \\
          output_dir/query/tree_species.treefile \\
          --rank species
""")


if __name__ == '__main__':
    # Run the demonstration
    demonstrate_all_ranks()
    
    # Show usage tips
    show_usage_tips()
