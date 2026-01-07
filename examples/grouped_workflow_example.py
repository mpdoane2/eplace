#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example script demonstrating the grouped ePLACE workflow for BLAST sequence comparison.

This script shows how to:
1. Read a query FASTA file
2. Run BLAST search against core_nt database
3. Filter results based on identity and coverage
4. Group sequences by a specified taxonomic rank (e.g., class, order)
5. Create one FASTA file per group containing all queries and unique references
6. Trim sequences to aligned regions
7. Build MAFFT alignments and IQTree phylogenies for each group

This workflow differs from blast_workflow_example.py in that it groups multiple
queries together based on their taxonomic classification (group_rank) and processes
them as a unit, rather than processing each query independently.
"""

import os
import sys
import argparse
import logging
from pathlib import Path

from eplace_lib.blast_analysis import run_blast_search, FastaReader
from eplace_lib.taxonomy import process_blast_results_for_taxonomy, rewrite_blast_hits, generate_classification_summary
from eplace_lib.alignment import check_alignment_consistency, group_hits_by_group_rank
from eplace_lib.alignment import create_grouped_fasta_with_queries, process_grouped_alignment_and_tree
from eplace_lib.alignment import process_grouped_alignment_and_tree_parallel, IQTreeBuilder

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Main function to run the grouped workflow."""
    parser = argparse.ArgumentParser(
        description='Run BLAST search and group sequences by taxonomic rank for joint analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters (group by class)
  python grouped_workflow_example.py query.fasta output_dir

  # Group by order instead of class
  python grouped_workflow_example.py query.fasta output_dir --group-rank order

  # Specify different rank for representatives and custom thresholds
  python grouped_workflow_example.py query.fasta output_dir --rank genus --group-rank family --min-identity 95

  # Use custom BLAST database location
  python grouped_workflow_example.py query.fasta output_dir --blastdb-path /path/to/blastdb

Notes:
  - This script requires BLAST+ tools (blastn and blastdbcmd) to be installed
  - The core_nt database must be downloaded first (see download_ncbi_example.py)
  - Default filtering: 90% identity over 80% query coverage
  - Outputs one set of files per taxonomic group (class, order, etc.)
  - Each group includes all queries that match to that taxonomic group
        """
    )
    
    parser.add_argument(
        'query_fasta',
        type=Path,
        help='Path to query FASTA file containing sequences to search'
    )
    
    parser.add_argument(
        'output_dir',
        type=Path,
        help='Directory for output files (will be created if it does not exist)'
    )
    
    parser.add_argument(
        '--rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for representative selection (default: genus)'
    )

    parser.add_argument(
        '--group-rank',
        type=str,
        default='class',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank to use for grouping sequences across queries (default: class)'
    )

    parser.add_argument(
        '--tree-label-rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank to use for tree labeling (default: genus)'
    )
    
    parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum percent identity for BLAST hits (default: 90.0)'
    )
    
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=80.0,
        help='Minimum query coverage percentage for BLAST hits (default: 80.0)'
    )
    
    parser.add_argument(
        '--database',
        type=str,
        default='core_nt',
        help='BLAST database name (default: core_nt)'
    )
    
    parser.add_argument(
        '--blastdb-path',
        type=Path,
        default=None,
        help='Path to BLAST database directory (default: $BLASTDB or ~/blastdb)'
    )
    
    parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='Number of threads for BLAST search and alignments (default: 1)'
    )

    parser.add_argument(
        '--overwrite-existing-blast',
        action='store_true',
        help='If the blast results already exist, redo the search and overwrite them'
    )
    
    parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip the alignment and tree building steps'
    )

    parser.add_argument(
        '--alignment-tolerance',
        type=int,
        default=50,
        help='Maximum allowed coordinate difference for alignment consistency check (default: 50)'
    )
    
    parser.add_argument(
        '--output-classification',
        type=Path,
        default=None,
        help='Path to output classification TSV file (default: <query_basename>_classification.tsv)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Display what would be done without actually running BLAST'
    )
    
    args = parser.parse_args()
    
    # Set default output classification file if not provided
    if args.output_classification is None:
        # Get the base name of the query fasta without extension
        base_name = args.query_fasta.stem
        # Remove common FASTA extensions if present
        for ext in ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']:
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
                break
        args.output_classification = args.output_dir / f"{base_name}_classification.tsv"
    if not args.output_classification.is_absolute():
        args.output_classification = args.output_dir / args.output_classification
    
    # Validate input file
    if not args.query_fasta.exists():
        logger.error(f"Query FASTA file not found: {args.query_fasta}")
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    skip_existing = not args.overwrite_existing_blast

    logger.info("=" * 60)
    logger.info("ePLACE Grouped BLAST Workflow")
    logger.info("=" * 60)
    logger.info(f"Query FASTA: {args.query_fasta}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Representative rank: {args.rank}")
    logger.info(f"Grouping rank: {args.group_rank}")
    logger.info(f"Tree labeling rank: {args.tree_label_rank}")
    logger.info(f"Classification output file: {args.output_classification}")
    logger.info(f"Min identity: {args.min_identity}%")
    logger.info(f"Min coverage: {args.min_coverage}%")
    logger.info(f"Overwrite: {args.overwrite_existing_blast} (skip_existing: {skip_existing})")
    logger.info(f"Database: {args.database}")
    logger.info(f"Threads: {args.num_threads}")
    
    if args.dry_run:
        logger.info("\nDRY RUN MODE - No BLAST search will be performed")
        logger.info("=" * 60)
        return 0
    
    logger.info("=" * 60)
    
    # Step 1: Read query sequences
    logger.info("\n[Step 1/7] Reading query sequences...")
    try:
        sequences = FastaReader.read_fasta(args.query_fasta)
        logger.info(f"Found {len(sequences)} query sequences")
        for seq_id in sequences.keys():
            logger.info(f"  - {seq_id} ({len(sequences[seq_id])} bp)")
    except Exception:
        logger.exception(f"Error reading FASTA file")
        sys.exit(1)
    
    # Step 2: Run BLAST search
    logger.info("\n[Step 2/7] Running BLAST search...")
    blast_output = args.output_dir / "blast_results.txt"
    
    try:
        success, filtered_hits = run_blast_search(
            query_fasta=args.query_fasta,
            output_file=blast_output,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            database=args.database,
            blastdb_path=args.blastdb_path,
            num_threads=args.num_threads,
            skip_existing=skip_existing
        )
        
        if not success:
            logger.error("BLAST search failed")
            sys.exit(1)
        
        logger.info(f"BLAST search completed successfully")
        logger.info(f"Found {len(filtered_hits)} hits after filtering")
        
    except Exception:
        logger.exception(f"Error during BLAST search")
        sys.exit(1)
    
    # Step 3: Process taxonomy information
    logger.info(f"\n[Step 3/7] Processing taxonomy information (rank: {args.rank})...")
    
    try:
        # This adds taxonomic information to the hits
        results = process_blast_results_for_taxonomy(
            blast_hits=filtered_hits,
            output_dir=args.output_dir,
            rank=args.rank,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        logger.info("\nExtracted representative sequences per query:")
        for query_id, output_fasta in results.items():
            if output_fasta:
                logger.info(f"  {query_id}: {output_fasta}")
            else:
                logger.warning(f"  {query_id}: Failed to extract sequences")
        
    except Exception:
        logger.exception(f"Error extracting sequences (grouped workflow)")
        sys.exit(1)

    logger.info("Rewriting the BLAST output file with the new annotations")
    try:
        rewrite_blast_hits(
            blast_hits=filtered_hits,
            output_file=args.output_dir / "blast_results_annotated.txt",
            header=True
        )
    except Exception as e:
        logger.exception(f"Error rewriting the blast hits: {e}")
        sys.exit(1)

    # Step 4: Check alignment consistency
    logger.info("\n[Step 4/7] Checking alignment consistency...")
    consistency = check_alignment_consistency(
        blast_hits=filtered_hits,
        tolerance=args.alignment_tolerance
    )
    
    inconsistent_count = sum(1 for is_consistent in consistency.values() if not is_consistent)
    if inconsistent_count > 0:
        logger.warning(
            f"Found {inconsistent_count} reference sequences with inconsistent alignments "
            f"(tolerance: {args.alignment_tolerance} bp)"
        )
    else:
        logger.info("All reference sequences have consistent alignments across queries")

    # Step 5: Group hits by group_rank
    logger.info(f"\n[Step 5/7] Grouping hits by {args.group_rank}...")
    grouped_hits = group_hits_by_group_rank(filtered_hits, args.group_rank)
    
    if not grouped_hits:
        logger.error("No groups found after grouping by rank")
        sys.exit(1)

    # Step 6: Create grouped FASTA files
    logger.info(f"\n[Step 6/7] Creating grouped FASTA files...")
    
    group_results = {}
    for group_tid, query_hits_map in grouped_hits.items():
        # Use the grouping key directly as the group name
        group_name = group_tid
        
        if not group_name:
            logger.warning(f"No group name found for group {group_tid}, skipping")
            continue
        
        logger.info(f"\nProcessing group: {group_name} ({group_tid})")
        logger.info(f"  Queries in this group: {', '.join(query_hits_map.keys())}")
        
        # Create group directory
        safe_group_name = group_name.replace(' ', '_').replace('/', '_').replace('|', '_')
        group_dir = args.output_dir / safe_group_name
        group_dir.mkdir(parents=True, exist_ok=True)
        
        # Create combined FASTA file
        combined_fasta = group_dir / f"{safe_group_name}_combined.fasta"
        success = create_grouped_fasta_with_queries(
            group_tid=group_tid,
            group_name=group_name,
            query_hits_map=query_hits_map,
            labeling_rank=args.rank,
            query_fasta=args.query_fasta,
            output_fasta=combined_fasta,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        if success:
            logger.info(f"  ✓ Created: {combined_fasta}")
            group_results[group_tid] = {
                'name': group_name,
                'dir': group_dir,
                'fasta': combined_fasta,
                'query_ids': list(query_hits_map.keys()),
                'hits': [hit for hits in query_hits_map.values() for hit in hits]
            }
        else:
            logger.error(f"  ✗ Failed to create grouped FASTA for {group_name}")

    # Step 7: Build alignments and trees for each group
    if not args.skip_alignment:
        logger.info("\n[Step 7/7] Building alignments and phylogenetic trees for each group...")
        
        # First, process all groups to do trimming and alignment
        # and start tree building in background
        tree_jobs = []
        group_job_info = {}
        
        for group_tid, group_info in group_results.items():
            group_name = group_info['name']
            group_dir = group_info['dir']
            query_ids = group_info['query_ids']
            blast_hits = group_info['hits']
            
            logger.info(f"\nProcessing group: {group_name} ({group_tid})")
            
            try:
                result = process_grouped_alignment_and_tree_parallel(
                    group_name=group_name,
                    group_dir=group_dir,
                    taxonomic_rank=args.tree_label_rank,
                    blast_hits=blast_hits,
                    query_ids=query_ids,
                    num_threads=args.num_threads,
                    background_tree=True  # Start tree building in background
                )
                
                # Log progress
                if result['alignment']:
                    logger.info(f"  ✓ Alignment: {result['alignment']}")
                if result['trimmed_fasta']:
                    logger.info(f"  ✓ Trimmed sequences: {result['trimmed_fasta']}")
                
                # Collect tree jobs for later waiting
                if result['tree_job']:
                    tree_jobs.append(result['tree_job'])
                    group_job_info[str(result['tree_file'])] = {
                        'group_name': group_name,
                        'tree_file': result['tree_file'],
                        'labeled_tree_path': result['labeled_tree_path'],
                        'blast_hits': result['blast_hits'],
                        'taxonomic_rank': result['taxonomic_rank']
                    }
                    logger.info(f"  ✓ Tree building started in background")
                    
            except Exception as e:
                logger.error(f"Error processing {group_name}: {e}")
                continue
        
        # Wait for all tree building jobs to complete
        if tree_jobs:
            logger.info(f"\n{'='*60}")
            logger.info(f"Waiting for {len(tree_jobs)} tree building jobs to complete...")
            logger.info(f"This may take a while depending on the size of the alignments.")
            logger.info(f"{'='*60}\n")
            
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            
            # Relabel trees that completed successfully
            logger.info("\nRelabeling trees with taxonomic names...")
            for tree_path, success in tree_results.items():
                if success and tree_path in group_job_info:
                    job_info = group_job_info[tree_path]
                    group_name = job_info['group_name']
                    tree_file = job_info['tree_file']
                    labeled_tree = job_info['labeled_tree_path']
                    blast_hits = job_info['blast_hits']
                    taxonomic_rank = job_info['taxonomic_rank']
                    
                    logger.info(f"  Processing {group_name}...")
                    if IQTreeBuilder.relabel_tree_with_taxonomy(
                        tree_file=tree_file,
                        blast_hits=blast_hits,
                        output_tree=labeled_tree,
                        taxonomic_rank=taxonomic_rank
                    ):
                        logger.info(f"    ✓ Labeled tree: {labeled_tree}")
                    else:
                        logger.warning(f"    ! Failed to relabel tree, but unlabeled tree is available: {tree_file}")
        
        logger.info(f"\nAlignment and tree building completed for {len(group_results)} groups")
    else:
        logger.info("\n[Step 7/7] Skipping alignment and tree building (--skip-alignment)")

    # Step 8: Generate classification summary TSV
    logger.info("\nGenerating classification summary TSV file...")
    try:
        success = generate_classification_summary(
            sequences=sequences,
            blast_hits=filtered_hits,
            output_file=args.output_classification,
            rank=args.rank,
            group_rank=args.group_rank,
            tree_label_rank=args.tree_label_rank
        )
        
        if success:
            logger.info(f"✓ Classification summary: {args.output_classification}")
        else:
            logger.warning("Failed to generate classification summary")
    except Exception as e:
        logger.error(f"Error generating classification summary: {e}")

    logger.info("\n" + "=" * 60)
    logger.info("Grouped workflow completed successfully!")
    logger.info(f"Processed {len(group_results)} taxonomic groups")
    logger.info("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
