#!/usr/bin/env python
"""
Example script demonstrating the complete ePLACE workflow for BLAST sequence comparison.

This script shows how to:
1. Read a query FASTA file
2. Run BLAST search against core_nt database
3. Filter results based on identity and coverage
4. Extract representative sequences per taxonomic rank
5. Save results to separate FASTA files
"""

import sys
import argparse
import logging
from pathlib import Path
from collections import defaultdict

from eplace_lib.blast_analysis import run_blast_search, FastaReader
from eplace_lib.taxonomy import process_blast_results_for_taxonomy, rewrite_blast_hits
from eplace_lib.alignment import process_query_alignment_and_tree

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Main function to run the complete workflow."""
    parser = argparse.ArgumentParser(
        description='Run BLAST search and extract representative sequences by taxonomic rank',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  python blast_workflow_example.py query.fasta output_dir

  # Specify taxonomic rank and custom thresholds
  python blast_workflow_example.py query.fasta output_dir --rank genus --min-identity 95 --min-coverage 85

  # Use custom BLAST database location
  python blast_workflow_example.py query.fasta output_dir --blastdb-path /path/to/blastdb

Notes:
  - This script requires BLAST+ tools (blastn and blastdbcmd) to be installed
  - The core_nt database must be downloaded first (see download_ncbi_example.py)
  - Default filtering: 90% identity over 80% query coverage
  - Outputs one FASTA file per query sequence in separate folders
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
        help='Taxonomic rank for representative selection (default: species)'
    )

    parser.add_argument(
        '--group-rank',
        type=str,
        default='class',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank to use to group the sequences across samples'
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
        help='Number of threads for BLAST search (default: 1)'
    )

    parser.add_argument(
        '--overwrite_existing_blast',
        action='store_true',
        help='If the blast results already exist, redo the search and overwrite them'
    )
    
    parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip the alignment and tree building steps'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Display what would be done without actually running BLAST'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not args.query_fasta.exists():
        logger.error(f"Query FASTA file not found: {args.query_fasta}")
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    skip_existing = not args.overwrite_existing_blast

    logger.info("=" * 60)
    logger.info("ePLACE BLAST Workflow")
    logger.info("=" * 60)
    logger.info(f"Query FASTA: {args.query_fasta}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Taxonomic rank: {args.rank}")
    logger.info(f"Taxonomic rank for grouping: {args.group_rank}")
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
    logger.info("\n[Step 1/5] Reading query sequences...")
    try:
        sequences = FastaReader.read_fasta(args.query_fasta)
        logger.info(f"Found {len(sequences)} query sequences")
        for seq_id in sequences.keys():
            logger.info(f"  - {seq_id} ({len(sequences[seq_id])} bp)")
    except Exception:
        logger.exception(f"Error reading FASTA file")
        sys.exit(1)
    
    # Step 2: Run BLAST search
    logger.info("\n[Step 2/5] Running BLAST search...")
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
    
    # Step 3: Group hits by query and display summary
    logger.info("\n[Step 3/5] Analyzing BLAST results...")
    hits_by_query = defaultdict(int)
    for hit in filtered_hits:
        hits_by_query[hit.query_id] += 1
    
    for query_id, count in hits_by_query.items():
        logger.info(f"  {query_id}: {count} hits")
    
    # Step 4: Extract representative sequences
    logger.info(f"\n[Step 4/5] Extracting representative sequences (rank: {args.rank})...")
    
    try:
        results = process_blast_results_for_taxonomy(
            blast_hits=filtered_hits,
            output_dir=args.output_dir,
            rank=args.rank,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        logger.info("\nResults:")
        for query_id, output_fasta in results.items():
            if output_fasta:
                logger.info(f"  {query_id}: {output_fasta}")
            else:
                logger.warning(f"  {query_id}: Failed to extract sequences")
        
    except Exception:
        logger.exception(f"Error extracting sequences (blast workflow)")
        sys.exit(1)

    logger.info("Rewriting the blast output file with the new annotations")
    try:
        rewrite_blast_hits(
            blast_hits=filtered_hits,
            output_file=args.output_dir / "blast_results_annotated.txt",
            header=True
        )
    except Exception:
        logger.exception(f"Error rewriting the blast hits")
        sys.exit(1)

    # Step 5: Align sequences and build trees (if not skipped)
    if not args.skip_alignment:
        logger.info("\n[Step 5/5] Aligning sequences and building phylogenetic trees...")
        
        # Group hits by query for processing
        hits_by_query_map = defaultdict(list)
        for hit in filtered_hits:
            hits_by_query_map[hit.query_id].append(hit)
        
        alignment_results = {}
        for query_id, query_hits in hits_by_query_map.items():
            logger.info(f"\nProcessing alignment and tree for query: {query_id}")
            
            # Get the query directory
            safe_query_id = query_id.replace('|', '_').replace('/', '_')
            query_dir = args.output_dir / safe_query_id
            
            if not query_dir.exists():
                logger.warning(f"Query directory not found: {query_dir}")
                continue
            
            # Process alignment and tree
            try:
                result = process_query_alignment_and_tree(
                    query_id=query_id,
                    query_dir=query_dir,
                    blast_hits=query_hits,
                    taxonomic_rank=args.rank,
                    query_fasta=args.query_fasta,
                    num_threads=args.num_threads
                )
                alignment_results[query_id] = result
                
                # Log results
                if result['labeled_tree']:
                    logger.info(f"  ✓ Labeled tree: {result['labeled_tree']}")
                elif result['tree']:
                    logger.info(f"  ✓ Tree: {result['tree']}")
                if result['alignment']:
                    logger.info(f"  ✓ Alignment: {result['alignment']}")
                if result['trimmed_fasta']:
                    logger.info(f"  ✓ Trimmed sequences: {result['trimmed_fasta']}")
                    
            except Exception as e:
                logger.error(f"Error processing {query_id}: {e}")
                continue
        
        logger.info(f"\nAlignment and tree building completed for {len(alignment_results)} queries")
    else:
        logger.info("\n[Step 5/5] Skipping alignment and tree building (--skip-alignment)")

    logger.info("\n" + "=" * 60)
    logger.info("Workflow completed successfully!")
    logger.info("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
